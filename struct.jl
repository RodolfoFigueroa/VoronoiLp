module DStruct
import Base.show
using Plots, Statistics
export DCEL, Vertex, Edge, Face, Handler,

voronoihelper, fixids!, plotdcel, checkdcel, mergevoronoi, reset, voronoi

#---
TOTAL_EDGES = 0
TOTAL_FACES = 0
TOTAL_VERTICES = 0

function reset()
    global TOTAL_EDGES = 0
    global TOTAL_FACES = 0
    global TOTAL_VERTICES = 0
    TOTAL_VERTICES = TOTAL_FACES = TOTAL_EDGES = 0
    return
end
#--

#Struct
mutable struct Vertex
    id::Int
    pos::Array
    edge
    original::Bool
    dead::Bool
    ghost::Bool
end
Vertex(id)::Vertex = Vertex(id, [NaN, NaN], nothing, true, false, false)
Vertex(id, pos)::Vertex = Vertex(id, pos, nothing, true, false, false)
Vertex(id, pos, edge)::Vertex = Vertex(id, pos, edge, true, false, false)
Vertex(id, pos, edge, original)::Vertex = Vertex(id, pos, edge, original, false, false)

mutable struct Face
    id::Int
    edge
    site::Union{Array,Nothing}
end
Face(id)::Face = Face(id, nothing, nothing)
Face(id, edge)::Face = Face(id, edge, nothing)
Face(id, edge, site)::Face = Face(id, edge, site)

const FI = Face(0)


mutable struct Edge
    id::Int
    orig::Union{Vertex,Nothing}
    dest::Union{Vertex,Nothing}
    cwo::Union{Edge,Nothing}
    ccwo::Union{Edge,Nothing}
    cwd::Union{Edge,Nothing}
    ccwd::Union{Edge,Nothing}
    fr::Union{Face,Nothing}
    fl::Union{Face,Nothing}
    dead::Bool
end
Edge(id)::Edge = Edge(id, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, false)
Edge(id, orig, dest)::Edge = Edge(id, orig, dest, nothing, nothing, nothing, nothing, nothing, nothing, false)


mutable struct DCEL
    edgelist::Array{Edge}
    vertexlist::Array{Vertex}
    facelist::Array{Face}
end
DCEL() = DCEL([], [], [])


mutable struct Handler
    left_face::Union{Face,Nothing}
    right_face::Union{Face,Nothing}
    left_vertex::Union{Vertex,Nothing}
    right_vertex::Union{Vertex,Nothing}
    current_vertex::Union{Vertex,Nothing}
    current_joint::Union{Edge,Nothing}
    ignore::Tuple
end
Handler(lf::Face, rf::Face, lv::Vertex, rv::Vertex)::Handler = Handler(lf, rf, lv, rv, nothing, nothing, ())



#---
#Printing
"""
    printnothing(x)

Print "e", "f" or "v" if `x` is an edge, vertex or face, respectively, alongside its id.
If `x` is nothing, print "nothing".
"""
function printnothing(e::Edge)::String
    return isnothing(e) ? "nothing" : "e$(e.id)"
end

function printnothing(e::Face)::String
    return isnothing(e) ? "nothing" : "f$(e.id)"
end

function printnothing(e::Vertex)::String
    return isnothing(e) ? "nothing" : "v$(e.id)"
end

function printnothing(e::Nothing)::String
    return "nothing"
end

function show(io::IO, v::Vertex)::Nothing
	    print(io, "v$(v.id): e=", printnothing(v.edge))
    print(io, " pos=$(v.pos) orig=$(v.original) ghost=$(v.ghost)")
    return
end


function show(io::IO, e::Edge)::Nothing
    print(io, "e$(e.id): orig=v$(e.orig.id) dest=v$(e.dest.id)")
    print(io, " cwo=", printnothing(e.cwo))
    print(io, " ccwo=", printnothing(e.ccwo))
    print(io, " cwd=", printnothing(e.cwd))
    print(io, " ccwd=", printnothing(e.ccwd))
    print(io, " fr=", printnothing(e.fr))
    print(io, " fl=", printnothing(e.fl))
    return
end


function show(io::IO, f::Face)::Nothing
    print(io, "f$(f.id): e=", printnothing(f.edge))
    print(io, " site=", isnothing(f.site) ? "nothing" : f.site)
    return
end


function show(io::IO, D::DCEL)::Nothing
    for v in D.vertexlist
        show(io, v)
        println(io)
    end
    println(io)
    for e in D.edgelist
        show(io, e)
        println(io)
    end
    println(io)
    for f in D.facelist
        show(io, f)
        println(io)
    end
    return
end



#---
#Geometry
"""
    cross2d(u, v)

Compute the 2-dimensional cross product u₁v₂-u₂v₁.
"""
@inline function cross2d(u::T, v::T) where T<:Union{Tuple, Array}
    return u[1]*v[2] - u[2]*v[1]
end


"""
    dot(u, v)

Compute the 2-dimensional dot product u₁v₁+u₂v₂.
"""
@inline function dot(u::Array, v::Array)::Float64
    return u[1]*v[1] + u[2]*v[2]
end


"""
    norm(x)


"""
@inline function norm(a::Array)::Float64
    return hypot(a...)
end


@inline function distance(p::Array, q::Array)::Float64
	return norm(p .-q)
end


@inline function distance(p::Vertex, q::Vertex)::Float64
	return distance(p.pos, q.pos)
end


"""
    cosatan(y, x)

Compute `cos(atan(y, x))`.

"""
@inline function cosatan(y::Number, x::Number)::Float64
    return x/hypot(x,y)
end


"""
    sinatan(y, x)

Compute `sin(atan(y, x))`.
"""
function sinatan(y::Number, x::Number)::Float64
    if x==0
        if y>0
            return 1
        elseif y<0
            return -1
        else
            return 0
        end
    else
        return y/hypot(x,y)
    end
    return
end


"""
    pointccw(x)

Tell whether the points in `x` are sorted counterclockwise.
"""
function pointccw(array::Array)::Bool
    sum = 0
    for i in 2:length(array)
        sum += (array[i][1] - array[i-1][1])*(array[i][2] + array[i-1][2])
    end
    sum += (array[1][1] - array[end][1])*(array[1][2] + array[end][2])
    return sum <= 0
end


function angleccw(a::Number, b::Number, c::Number)::Bool
    return sin(a-b) + sin(b-c) + sin(c-a)  <= 0
end


"""
    circlethreepoints(x, y, z)

Return the center of the circle passing through the points `x`, `y` and `z`.
"""
function circlethreepoints(p::Array, q::Array, r::Array)::Array
    x1,y1 = p
    x2,y2 = q
    x3,y3 = r
    den = 2*(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))
    x = (x1^2+y1^2)*(y3-y2)+(x2^2+y2^2)*(y1-y3)+(x3^2+y3^2)*(y2-y1)
    y = (x1^2+y1^2)*(x2-x3)+(x2^2+y2^2)*(x3-x1)+(x3^2+y3^2)*(x1-x2)
    return [x/den,y/den]
end


function circlethreepoints(p::Face, q::Face, r::Face)::Array
	return circlethreepoints(p.site, q.site, r.site)
end



#--- Distances
function l2(x::Array, y::Array)::Float64
	return hypot(x,y)
end


#---
"""
    settopology(e)

"""
function settopology!(e::Edge; orig::Union{Vertex,Nothing}=nothing,
	dest::Union{Vertex,Nothing}=nothing, cwo::Union{Edge,Nothing}=nothing,
    ccwo::Union{Edge,Nothing}=nothing, cwd::Union{Edge,Nothing}=nothing,
    ccwd::Union{Edge,Nothing}=nothing, fr::Union{Face,Nothing}=nothing,
    fl::Union{Face,Nothing}=nothing)::Nothing
	if !isnothing(orig) e.orig = orig end
	if !isnothing(dest) e.dest = dest end
	if !isnothing(cwo) e.cwo = cwo end
    if !isnothing(ccwo) e.ccwo = ccwo end
    if !isnothing(cwd) e.cwd = cwd end
    if !isnothing(ccwd) e.ccwd = ccwd end
    if !isnothing(fr) e.fr = fr end
    if !isnothing(fl) e.fl = fl end
    return
end
settopology!(e::Edge)::Nothing = settopology!(e, cwo=e, ccwo=e, cwd=e, ccwd=e)



#---
"""
    createvertex!(D, pos, orig=true)

Create a vertex with position `pos` and original flag `orig` and add it to the DCEL `D`.
"""
function createvertex!(D::DCEL, pos::Array, original::Bool=true)::Vertex
    global TOTAL_VERTICES
    TOTAL_VERTICES += 1
    new_vertex = Vertex(TOTAL_VERTICES, pos, nothing, original)
    push!(D.vertexlist, new_vertex)
    return new_vertex
end


function createghostvertex!(D::DCEL, pos::Array)::Vertex
    v = createvertex!(D, pos)
    v.ghost = true
    return v
end


"""
    createedge!(D, p, q)

Create an edge with origin `p` and destination `q` and add it to the DCEL `D`.
"""
function createedge!(D::DCEL, u::Vertex, v::Vertex)::Edge
    global TOTAL_EDGES
    TOTAL_EDGES += 1
    new_edge = Edge(TOTAL_EDGES, u, v)
    settopology!(new_edge)
    push!(D.edgelist, new_edge)
    return new_edge
end


"""
    createface!(D, e, s)

Create a face with edge `e` and site `s` and add it to the DCEL `D`.
"""
function createface!(D::DCEL, e::Edge, s::Array)
    global TOTAL_FACES
    TOTAL_FACES += 1
    f = Face(TOTAL_FACES, e, s)
    push!(D.facelist, f)
    return f
end


"""
    squeezeedge(v, e; <keyword_arguments>)

Update all the edges incident to `v` so that they are properly sorted around edge `e`.

`e` must have `v` as one of its endpoints. Checking if the edges are sorted correctly
can be skipped if a previous/next edge is provided via keyword arguments.

# Arguments
- `v::Vertex` : The vertex around which to sort.
- `e::Edge` : The edge to consider when sorting.
- `update_faces::Bool=true`: Whether to automatically update `e.fl` and `e.fr` with the proper faces.
- `previous::Edge`: The edge that will be clockwise to `e` around `v`.
- `next::Edge`: The edge that will be counterclockwise to `e` around `v`.
"""
function squeezeedge!(v::Vertex, e::Edge; update_faces::Bool=true,
	previous::Union{Edge,Nothing}=nothing, next::Union{Edge,Nothing}=nothing)::Nothing
	# @assert v in endpoints(e)
    if isnothing(v.edge)
        return
    end
    previous_set = next_set = true
    if isnothing(previous) && isnothing(next)
        previous_set = next_set = false
        previous = v.edge
    elseif isnothing(previous)
        previous_set = false
        previous = cw(next, v)
    elseif isnothing(next)
        next_set = false
    end
    while true
        if !next_set
            next = ccw(previous, v)
        end
        if previous_set || next_set || sorted_ccw(previous, e, next)
			fr = fl = nothing
			if update_faces
	            fr = v==previous.orig ? previous.fl : previous.fr
	            fl = v==next.orig ? next.fr : next.fl
			end
            if v == e.orig
                settopology!(e, cwo=previous, ccwo=next, fr=fr, fl=fl)
            else
                settopology!(e, cwd=previous, ccwd=next, fr=fl, fl=fr)
            end
            cwset!(next, v, e)
            ccwset!(previous, v, e)
            return
        end
        previous = next
    end
    return
end


"""
    joinvertices!(D, p, q; <keyword_arguments>)

Join the vertices `p` and `q`, add the resulting edge and faces to `D` and update the topology.

# Arguments
- `D::DCEL`: The DCEL that contains `p` and `q`.
- `p::Vertex`: The origin vertex to join.
- `q::Vertex`: The destination vertex to join.
- `previous_orig::Edge`: Edge that will be clockwise to the new edge around `p`.
- `next_orig::Edge`: Edge that will be counterclockwise to the new edge around `p`.
- `previous_dest::Edge`: Edge that will be clockwise to the new edge around `q`.
- `next_dest::Edge`: Edge that will be counterclockwise to the new edge around `q`.
- `split_face::Bool`: Whether to split the face that contains the new edge in two new faces.
- `update_edges::Bool`: Whether to update the vertices so that their assigned edge is the newly created edge.
- `skip_checks::Bool`: If set to true, sets all the edge fields of the new edge to itself and returns.
"""
function joinvertices!(D::DCEL, u::Vertex, v::Vertex; previous_orig::Union{Edge,Nothing}=nothing, 
    previous_dest::Union{Edge,Nothing}=nothing, next_orig::Union{Edge,Nothing}=nothing, 
    next_dest::Union{Edge,Nothing}=nothing, split_face::Bool=true, update_edges::Bool=true,
    skip_checks::Bool=false)
        
    new_edge = createedge!(D, u, v)
    if skip_checks
        settopology!(new_edge)
        u.edge = v.edge = new_edge
        return new_edge
    end
    squeezeedge!(u, new_edge, previous=previous_orig, next=next_orig)
    squeezeedge!(v, new_edge, previous=previous_dest, next=next_dest)
    disconnected = isnothing(u.edge) || isnothing(v.edge)
	if update_edges
    	u.edge = v.edge = new_edge
	end
    if disconnected
        return new_edge
    end

    if split_face
        global TOTAL_FACES
        f = new_edge.fr
        f1 = Face(TOTAL_FACES+1, new_edge)
        current_vertex = new_edge.orig
        current_edge = ccw(new_edge, current_vertex)
        while current_edge != new_edge
            if current_vertex == current_edge.orig
                current_edge.fr = f1
                current_vertex = current_edge.dest
            else
                current_edge.fl = f1
                current_vertex = current_edge.orig
            end
            current_edge = ccw(current_edge, current_vertex)
        end

        f2 = Face(TOTAL_FACES+2, new_edge)
        current_vertex = new_edge.dest
        current_edge = ccw(new_edge, current_vertex)
        while current_edge != new_edge
            if current_vertex == current_edge.orig
                current_edge.fr = f2
                current_vertex = current_edge.dest
            else
                current_edge.fl = f2
                current_vertex = current_edge.orig
            end
            current_edge = ccw(current_edge, current_vertex)
        end

        new_edge.fl = f1
        new_edge.fr = f2

        if !isnothing(f) && !isnothing(f.site)
            leftofedge(new_edge, f.site) ? new_edge.fl.site = f.site : new_edge.fr.site = f.site
        end
        
        push!(D.facelist, f1, f2)
        TOTAL_FACES += 2
    end
    return new_edge
end


"""
    createfloatingedge(D)

Create an edge with null origin and destination and add it to the DCEL D.
"""
function createfloatingedge!(D::DCEL)::Edge
    e = joinvertices!(D, Vertex(0), Vertex(0), skip_checks=true)
	return e
end


"""
    createdummyvertex!(D, u, angle)

Create a dummy vertex with an angle `angle` with respect to the vertex `u` and add
it to the DCEL D.

"""
function createdummyvertex!(D::DCEL, u::Vertex, angle::Number)::Vertex
    vector = u.pos .+ [cos(angle),sin(angle)]
    v = createvertex!(D, vector, false)
    return v
end


"""
    fixids!(D)

Fix the ids of the vertices, edges and faces of D so that they are sequential.
"""
function fixids!(D::DCEL)::Nothing
    for i in 1:length(D.vertexlist)
        D.vertexlist[i].id = i
    end
    for i in 1:length(D.edgelist)
        D.edgelist[i].id = i
    end
    for i in 1:length(D.facelist)
        D.facelist[i].id = i
    end
    return
end


"""
    averagevertexpositions(D)

Return the average of all the vertex positions of D.
"""
function averagevertexpositions(D::DCEL)
	pos = getfield.(D.vertexlist, :pos)
	mean_pos = mean(pos)
	dist = [abs.(p - mean_pos) for p in pos]
	mean_dist = mean(dist)
	return mean_pos, mean_dist
end


"""
    plotdcel(D; <keyword_arguments>)

Plot the DCEL D. 

# Arguments
- `faces::Bool`: 

Additionaly, 
"""
function plotdcel(D::DCEL; dead_edges::Bool=false, dead_vertices::Bool=false, 
                bounds::Bool=true, sites::Bool=false, normalize::Bool=true, 
                scale::Number=1, font_size::Int=12, face_labels::Bool=true,
                labels::Bool=true, line::Symbol=:arrow, show_all::Bool=false,
                safe_mode::Bool=false, kwargs...)::Plots.Plot

    if show_all
        p = plot(leg=false; kwargs...)
    else
        avg_pos, avg_dist = averagevertexpositions(D)
        xlims = avg_pos[1] .+ [-1, 1] .*avg_dist
        ylims = avg_pos[2] .+ [-1, 1] .*avg_dist
        p = plot(leg=false, xlims=xlims, ylims=ylims; kwargs...)
    end

	edge_pos_x = Array{Float64}(undef,0)
	edge_pos_y = Array{Float64}(undef,0)
	edge_color = []
    edge_text = []
    midpoint_x = []
    midpoint_y = []
    for e in D.edgelist
		start = finish = nothing
		if isstrut(e)
			start = e.orig.pos
			finish = e.dest.pos + (normalize ? normalizedummy(e.dest)*(scale-1) : [0,0])
		elseif isframe(e)
			start = e.orig.pos + (normalize ? normalizedummy(e.orig)*(scale-1) : [0,0])
			finish = e.dest.pos + (normalize ? normalizedummy(e.dest)*(scale-1) : [0,0])
		else
			start = e.orig.pos
			finish = e.dest.pos
		end
        bound = isframe(e)
        color = e.dead ? :gray : (bound ? :red : :black)
        ave = (start + finish)/2
        if (dead_edges || !e.dead) && (bounds || !bound)
            if !isnan(start[1]) && !isnan(finish[1])
                push!(edge_pos_x, start[1], finish[1], NaN)
                push!(edge_pos_y, start[2], finish[2], NaN)
                push!(midpoint_x, (start[1]+finish[1])/2)
                push!(midpoint_y, (start[2]+finish[2])/2)
                push!(edge_color, color)
                if labels
                    push!(edge_text, (ave[1], ave[2], "e$(e.id)", Plots.font("Sans",font_size)))
                end
            end
        end
    end
    plot!(p, edge_pos_x, edge_pos_y,color=edge_color, line=line, annotations=edge_text)

	vertex_pos_x = Array{Float64}(undef, 0)
	vertex_pos_y = Array{Float64}(undef, 0)
	vertex_col = Array{Symbol}(undef, 0)
	vertex_text = Array{Any}(undef, 0)
	for v in D.vertexlist
		# color = nothing
		if v.original
			color = :black
			pos = v.pos
		else
			color = :red
			pos = v.pos + (normalize ? normalizedummy(v)*(scale-1) : [0,0])
		end
		if v.dead
			color = :grey
		end

        if (dead_vertices || !v.dead) && (bounds || v.original)
            if !isnan(pos[1])
                push!(vertex_pos_x, pos[1])
                push!(vertex_pos_y, pos[2])
                push!(vertex_col, color)
                if labels
                    push!(vertex_text, (pos[1], pos[2], "v$(v.id)", Plots.font("Sans",font_size)))
                end
            end
		end
	end
	scatter!(p, vertex_pos_x, vertex_pos_y, color=vertex_col, annotations=vertex_text)

    face_text = []
    site_pos_x = []
    site_pos_y = []
    for f in D.facelist
        if face_labels
            vertices = faceverticesccw(f)
            x, y = mean(getfield.(vertices, :pos))
            push!(face_text, (x, y, "f$(f.id)", Plots.font("Sans",font_size)))
        end
        if sites && !isnothing(f.site)
            push!(site_pos_x, f.site[1])
            push!(site_pos_y, f.site[2])
        end
    end
    plot!(p, annotations=face_text)
    scatter!(p, site_pos_x, site_pos_y, color=:blue)
    return p
end


function checkedge(e::Edge)::Nothing
	@assert cw(e.ccwo, e.orig) == e "ccwo->cwo $(e)"
	@assert ccw(e.cwo, e.orig) == e "cwo->ccwo $(e)"
	@assert cw(e.ccwd, e.dest) == e "ccwd->cwd $(e)"
	@assert ccw(e.cwd, e.dest) == e "cwd->ccwd $(e)"
	@assert e.orig != e.dest
	if isstrut(e)
		@assert !e.dest.original
	end
	return
end


function writenothing(io, s::String)::Nothing
	if !isnothing(io)
		write(io, s)
	end
	return
end



function checkdcel(D::DCEL; io=stdout)::Nothing
    for v in D.vertexlist
		#writenothing(io, "Checking vertex: $v\n")
        @assert v in endpoints(v.edge)
    end
    for e in D.edgelist
		#writenothing(io, "Checking edge: $e\n")
        checkedge(e)
    end
    for f in D.facelist
		#writenothing(io, "Checking face: $f\n")
        edgesccw = faceedges(f, :ccw)
        edgescw = faceedges(f, :cw)
        @assert isempty(setdiff(edgescw, edgesccw))
        for e in edgesccw
            @assert f==e.fr || f==e.fl "$(f):$(e)"
        end
        @assert isempty(setdiff(faceverticescw(f), faceverticescw(f)))
    end
    return
end


function cw(e::Edge, v::Vertex)::Edge
	if v == e.orig
		return e.cwo
	elseif v == e.dest
		return e.cwd
	else #sanity
		throw("Looking for $v in e$(e.id) with orig=$(e.orig) and dest=$(e.dest)")
	end
end


function ccw(e::Edge, v::Vertex)::Edge
	if v == e.orig
		return e.ccwo
	elseif v == e.dest
		return e.ccwd
	else #sanity
		throw("Looking for $v in e$(e.id) with orig=$(e.orig) and dest=$(e.dest)")
	end
end


function nextpivot(e::Edge, v::Vertex, side::Symbol)::Edge
    # return getfield(Main, side)(e, v)
    if side == :cw
        return cw(e, v)
    else
        return ccw(e,v)
    end
    return
end


function cwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? settopology!(e, cwo=new_edge) : settopology!(e, cwd=new_edge)
    return
end


function ccwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? settopology!(e, ccwo=new_edge) : settopology!(e, ccwd=new_edge)
    return
end


@inline function endpoints(u::Edge)::Tuple
    return u.orig, u.dest
end


function edgeangle(u::Edge)::Number
    return atan(u.dest.pos[2]-u.orig.pos[2], u.dest.pos[1]-u.orig.pos[1])
end


function pivotangle(u::Edge, p::Vertex)::Number
    return edgeangle(u) + (p==u.orig ? 0 : pi)
end


function returnempty(s::Array)::Any
    return isempty(s) ? nothing : pop!(s)
end


# function commonvertex(u::Edge, edges...)::Union{Vertex,Nothing}
#     f = endpoints(u)
#     for e in endpoints.(edges)
#         intersect!(f, e)
#     end
#     return returnempty(f)
# end


function commonvertex(u::Edge, v::Edge)::Union{Vertex,Nothing}
	if u.orig == v.orig || u.orig == v.dest
		return u.orig
	elseif u.dest == v.orig || u.dest == v.dest
		return u.dest
	end
	return nothing
end


function commonvertex(u::Edge, v::Edge, w::Edge)::Union{Vertex,Nothing}
    a = u.orig
    b = u.dest
    c = v.orig
    d = v.dest
    e = w.orig
    f = w.dest
    if (a == c == e) || (a == c == f) || (a == d == e) || (a == d == f)
		return a
	elseif (b == c == e) || (b == c == f) || (b == d == e) || (b == d == f)
		return b
	end
	return
end

function uncommonvertices(u::Edge, v::Edge, p::Vertex)::Tuple
    return (u.orig==p ? u.dest : u.orig) , (v.orig==p ? v.dest : v.orig)
end


function uncommonvertices(u::Edge, v::Edge, w::Edge)::Tuple
    p = commonvertex(u, v, w)
    a = p == u.orig ? u.dest : u.orig
    b = p == v.orig ? v.dest : v.orig
    c = p == w.orig ? w.dest : w.orig
    return a, b, c, p
end


function sorted_ccw_1(u::Edge, v::Edge, w::Edge)::Bool
    pivot = commonvertex(u, v, w)
    a = pivotangle(u, pivot)
    b = pivotangle(v, pivot)
    c = pivotangle(w, pivot)
    return sin(a-b) + sin(b-c) + sin(c-a) <= 0
end


function sorted_ccw(u::Edge, v::Edge, w::Edge)::Bool
    a, b, c, p = getfield.(uncommonvertices(u, v, w), :pos)
    a -= p
    b -= p
    c -= p
    return cross2d(a, b) * sqrt(c[1]^2+c[2]^2) + cross2d(b, c) * sqrt(a[1]^2+a[2]^2) + cross2d(c, a) * sqrt(b[1]^2+b[2]^2) >= 0
end


function getstrut(v::Vertex)::Edge #TODO: add sanity
    @assert !v.original "v is not a dummy vertex"
    edge = v.edge
    while !isstrut(edge)
        edge = ccw(edge, v)
    end
    return edge
end


function isframe(e::Edge)::Bool
    # return all(.!endfield(e, :original))
	return !e.orig.original && !e.dest.original
end


function findframe(f::Face; dir::Symbol=:ccw)::Edge
	edge = f.edge
	while !isframe(edge)
		edge = nextedge(f, dir, edge)
	end
	return edge
end


function normalizedummy(v::Vertex)::Array
    @assert !v.original "v is not a dummy vertex"
    e = getstrut(v)
    u, v = e.orig.original ? (e.orig, e.dest) : (e.dest, e.orig)
    return v.pos - u.pos
end


vertexccw(vertices::Array) = pointccw(getfield.(vertices, :pos))


function isstrut(e::Edge)::Bool
    # a = endfield(e, :original)
    # return xor(a[1], a[2])
	return (e.orig.original && !e.dest.original) || (!e.orig.original && e.dest.original)
end


function endfield(e::Edge, f::Symbol)
    return getfield.(endpoints(e), f)
end


function cwface(f::Face, edge::Union{Edge,Nothing}=nothing)::Edge
    if isnothing(edge)
        return f.edge
    end
    # if edge.fr == edge.fl
    #     return edge.ccwo != edge ? edge.ccwo : edge.ccwd
    if f == edge.fr
        return edge.ccwd
    elseif f == edge.fl
        return edge.ccwo
    else #sanity
        throw("Looking for $(f) in e$(edge.id) with fl=$(edge.fl) and fr=$(edge.fr)")
    end
    return
end


function ccwface(f::Face, edge::Union{Edge,Nothing}=nothing)::Edge
    if isnothing(edge)
        return f.edge
    end
    # if edge.fr == edge.fl
    #     return edge.cwo != edge ? edge.cwo : edge.cwd
    if f == edge.fr
        return edge.cwo
    elseif f == edge.fl
        return edge.cwd
    else
        throw("Looking for $(f) in e$(edge.id) with fl=$(edge.fl) and fr=$(edge.fr)")
    end
    return
end


function nextedge(f::Face, dir::Symbol, start::Union{Edge,Nothing}=nothing)::Edge
	if dir == :ccw
		return ccwface(f, start)
	elseif dir == :cw
		return cwface(f, start)
	else
		throw("$dir is not a valid direction. Possible directions are clockwise (:cw) and counterclockwise (:ccw)")
	end
	return
end


function faceedges(f::Face, dir::Symbol, edge::Union{Edge,Nothing}=nothing)::Array
	if isnothing(edge)
		edge = f.edge
	end
    out = [edge]
    edge = nextedge(f, dir, edge)
	while edge != f.edge
		push!(out, edge)
		edge = nextedge(f, dir, edge)
	end
	return out
end


function faceverticescw(f::Face)::Array
    out = []
    edge = f.edge
    while true
        push!(out, f==edge.fr ? edge.orig : edge.dest)
        edge = ccwface(f, edge)
        if edge == f.edge
            break
        end
    end
    return out
end


function faceverticesccw(f::Face)::Array
    out = []
    edge = f.edge
    while true
        push!(out, f==edge.fr ? edge.orig : edge.dest)
        edge = cwface(f, edge)
        if edge == f.edge
            break
        end
    end
    return out
end


function facevertices(f::Face, dir::Symbol=:ccw)::Array
	if dir == :ccw
		return faceverticesccw(f)
	elseif dir == :cw
		return faceverticescw(f)
	else
		throw("$dir is not a valid direction. Possible directions are clockwise (:cw) and counterclockwise (:ccw)")
	end
	return
end


function vertexedges(v::Vertex, dir::Symbol=:ccw)::Array
    starter_edge = v.edge
    edge = nextpivot(starter_edge, v, dir)
	out = [starter_edge]
	while edge != starter_edge
        push!(out, edge)
        edge = nextpivot(edge, v, dir)
	end
	return out
end


function findfaceframe(f::Face, edge::Union{Edge,Nothing}=nothing; dir::Symbol=:ccw)::Union{Edge,Nothing}
    isnothing(edge) ? edge=f.edge : @assert f==edge.fr || f==edge.fl
    starting_edge = edge
    while !isframe(edge)
		edge = nextedge(f, dir, edge)
        if edge == starting_edge
            return nothing
        end
    end
    return edge
end


function splitedge!(D::DCEL, e::Edge, split::Vertex)::Tuple
    new_edge = Edge(e.id, split, e.dest)
    push!(D.edgelist, new_edge)

    cwd = e.cwd==e ? new_edge : e.cwd
    ccwd = e.ccwd==e ?  new_edge : e.ccwd
    settopology!(new_edge, cwo=e, ccwo=e, cwd=cwd, ccwd=ccwd, fr=e.fr, fl=e.fl)
    split.edge = new_edge
    cwset!(e.ccwd, e.dest, new_edge)
    ccwset!(e.cwd, e.dest, new_edge)
    if isstrut(e)
        x, y = e.dest.pos - e.orig.pos
        new_edge.dest.pos = split.pos + [cosatan(y,x), sinatan(y,x)]
    end

    e.dest.edge = new_edge
    e.dest = split
    e.cwd = e.ccwd = new_edge
    return e, new_edge
end


function leftofedge(edge::Edge, c::Array, distance::Function=l2)::Bool
	a = edge.fl.site
	b = edge.fr.site
	return distance(a, c) < distance(b, c)
end
leftofedge(e::Edge, v::Vertex)::Bool = leftofedge(e, v.pos)


function leftofline(a::Array, b::Array, c::Array)::Bool
    return ((b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1])) > 0
end


function midpoint(f1::Face, f2::Face)::Array
    return mean([f1.site, f2.site])
end


function perpangle(start::Array, finish::Array)::Float64
    return atan(finish[2]-start[2], finish[1]-start[1]) + pi/2
end
perpangle(f1::Face, f2::Face)::Float64 = perpangle(f1.site, f2.site)


function perpvector(start::Array, finish::Array)::Array
    v = finish - start
    return [-v[2], v[1]]/hypot(v[1], v[2])
end
perpvector(f1::Face, f2::Face)::Array = perpvector(f1.site, f2.site)


function findliveedge(v::Vertex, dir::Symbol=:ccw)
	edge = v.edge
	while edge.dead
		edge = nextpivot(edge, v, dir)
	end
	return edge
end


function unstickedge!(e::Edge, v::Vertex)::Nothing
    if v == e.orig
        ccwset!(e.cwo, e.orig, e.ccwo)
        cwset!(e.ccwo, e.orig, e.cwo)
    elseif v == e.dest
        ccwset!(e.cwd, e.dest, e.ccwd)
        cwset!(e.ccwd, e.dest, e.cwd)
    else
        throw("Given vertex is not in endpoints of edge")
    end
	if v.edge == e && ccw(e, v) != e
		v.edge = findliveedge(v)
	end
    return
end


function unstickedge!(e::Edge, s::Symbol)::Nothing
	return unstickedge!(e, getfield(e, s))
end


function bisectorangles(u::Array, v::Array, w::Array)::Tuple
    p = perpangle(v, w)
    q = perpangle(w, u)
    r = perpangle(u, v)
    return p, q, r
end


function voronoitwopoints(points::Array, p::Int=2)::DCEL # Handmade so it's faster
    u, v = points
    
    D = DCEL()
    # m = createvertex!(D, mean([u,v]))
    m = createghostvertex!(D, mean([u,v]))
    angle = perpangle(v, u)
    
    x = createdummyvertex!(D, m, angle)
    y = createdummyvertex!(D, m, angle+pi)

    e1 = joinvertices!(D, m, x, skip_checks=true)
    e2 = joinvertices!(D, m, y, skip_checks=true)
    e3 = joinvertices!(D, x, y, skip_checks=true)
    e4 = joinvertices!(D, y, x, skip_checks=true)
    
    f1 = createface!(D, e1, u)
    f2 = createface!(D, e2, v)

    squeezeedge!(m, e1, previous=e2, next=e2)

    squeezeedge!(x, e1, next=e4, previous=e3)
    squeezeedge!(y, e2, next=e3, previous=e4)

    squeezeedge!(x, e3, next=e1, previous=e4)
    squeezeedge!(y, e4, next=e2, previous=e3)

    settopology!(e1, fl=f2, fr=f1)
    settopology!(e2, fl=f1, fr=f2)
    settopology!(e3, fl=f2, fr=FI)
    settopology!(e4, fl=f1, fr=FI)

    return D
end


function voronoithreepoints(points::Array)::DCEL # Handmade so it's faster
    u, v, w = points
    u, v = pointccw(points) ? (v, u) : (u, v) # Swap the first two points if they are not sorted properly
    a, b, c = bisectorangles(u, v, w)

    D = DCEL()
    m = createvertex!(D, circlethreepoints(u, v, w))

    x = createdummyvertex!(D, m, a)
    y = createdummyvertex!(D, m, b)
    z = createdummyvertex!(D, m, c)

    e1 = joinvertices!(D, m, x, skip_checks=true)
    e2 = joinvertices!(D, m, y, skip_checks=true)
    e3 = joinvertices!(D, m, z, skip_checks=true)
    e4 = joinvertices!(D, x, y, skip_checks=true)
    e5 = joinvertices!(D, y, z, skip_checks=true)
    e6 = joinvertices!(D, z, x, skip_checks=true)

    f1 = createface!(D, e2, u)
    f2 = createface!(D, e3, v)
    f3 = createface!(D, e1, w)

    squeezeedge!(m, e1, previous=e2, next=e3)
    squeezeedge!(m, e2, previous=e3, next=e1)
    squeezeedge!(m, e3, previous=e1, next=e2)
    
    squeezeedge!(x, e1, previous=e6, next=e4)
    squeezeedge!(y, e2, previous=e4, next=e5)
    squeezeedge!(z, e3, previous=e5, next=e6)

    squeezeedge!(x, e4, previous=e1, next=e6)
    squeezeedge!(y, e5, previous=e2, next=e4)
    squeezeedge!(z, e6, previous=e3, next=e5)

    settopology!(e1, fl=f2, fr=f3)
    settopology!(e2, fl=f3, fr=f1)
    settopology!(e3, fl=f1, fr=f2)
    settopology!(e4, fl=FI, fr=f3)
    settopology!(e5, fl=FI, fr=f1)
    settopology!(e6, fl=FI, fr=f2)

    return D
end


function frameface(e::Edge)::Face
    if e.fr == FI
        return e.fl
    elseif e.fl == FI
        return e.fr
    else #sanity
        throw("")
    end
    return
end


function framesite(e::Edge)::Array
    return frameface(e).site
end


# function findside1(D::DCEL, max::Bool)
#     frames = findall(x-> isframe(x), D.edgelist)
#     edges = D.edgelist[frames]
#     f = max ? argmax : argmin
# 	return edges[f(hcat(getfield.(frameface.(edges), :site)...)[1,:])]
# end


function findside(D::DCEL, max::Bool) #TOOD: Is this faster than above?
    frame = D.edgelist[findfirst(x -> isframe(x), D.edgelist)]
    previous_frame = max ? cwface(FI, frame) : ccwface(FI, frame)
    next_frame = max ? ccwface(FI, frame) : cwface(FI, frame)
    
    site = framesite(frame)[1]
    previous_site = framesite(previous_frame)[1]
    next_site = framesite(next_frame)[1]
    if max
        while site < previous_site || site < next_site
            previous_frame, frame = frame, next_frame
            next_frame = ccwface(FI, frame)
            previous_site, site = site, next_site
            next_site = framesite(next_frame)[1]
        end
    else
        while site > previous_site || site > next_site
            previous_frame, frame = frame, next_frame
            next_frame = cwface(FI, frame)
            previous_site, site = site, next_site
            next_site = framesite(next_frame)[1]
        end
    end
    return frame
end


function findsupport(left::DCEL, right::DCEL, top::Bool)::Tuple #TODO: Optimize
    global FI

	edge_l = findside(left, true)
    edge_r = findside(right, false)
    f_l = top ? cwface : ccwface
    f_r = top ? ccwface : cwface
    next_l = f_l(FI, edge_l)
    next_r = f_r(FI, edge_r)
    cond = top ? !leftofline : leftofline
    while cond(framesite(edge_r), framesite(edge_l), framesite(next_l)) || cond(framesite(edge_r), framesite(edge_l), framesite(next_r))
        framesite(edge_r)
        framesite(edge_l)
        framesite(next_l)
        framesite(next_r)    
        while cond(framesite(edge_r), framesite(edge_l), framesite(next_l))
            edge_l = next_l
            next_l = f_l(FI, edge_l)
        end
        while cond(framesite(edge_r), framesite(edge_l), framesite(next_r))
            edge_r = next_r
            next_r = f_r(FI, edge_r)
        end
    end
    return frameface(edge_l), frameface(edge_r)
end


function findextrema(left::DCEL, right::DCEL)::Tuple
    return (findsupport(left, right, true)..., findsupport(left, right, false)...)
end


function bisector(f1::Face, f2::Face)::DCEL
    u = f1.site
    v = f2.site
    D = DCEL()
    m = createvertex!(D, mean([u,v]))
    angle = atan(v[2]-u[2], v[1]-u[1]) + pi/2
    addray!(D, m, angle)
    addray!(D, m, angle+pi)
    return D
end


function edgerayintersection(edge::Edge, q::Array, vector::Array, infinite::Bool=false; io=stdout)::Array
    p = edge.orig.pos
    r = edge.dest.pos - p
    num = q - p
    den = cross2d(r, vector)
    t = cross2d(num, vector)/den
    if (cross2d(num, r)/den>=0 || infinite) && (isstrut(edge) ? 0<=t : 0<=t<=1)
        return p + t*r
    else
        return [NaN, NaN]
    end
    return
end


function facerayintersection(f::Face, start::Array, vector::Array, edge::Edge;  
	dir::Symbol=:ccw, infinite::Bool=false, ignore::Tuple, io=stdout)::Tuple
	start_vertices = endpoints(edge)
	start_edge = edge
	counter = 0
    #writenothing(io, "CHECKING INTERSECTION OF FACE $(f) WITH IGNORE $(getfield.(ignore, :id))\n")
    while true
        #writenothing(io, "CHECKING EDGE $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        #writenothing(io, "LEFT FACE: $(edge.fl), RIGHT FACE: $(edge.fr)\n")
        if !edge.dead && !isframe(edge) && !(edge in ignore)
            inter = edgerayintersection(edge, start, vector, infinite, io=io)
			#writenothing(io, "INTERSECTION: $inter\n")
            if !isnan(inter[1])
                #writenothing(io, "INTERSECTED $(edge.id) ($(edge.orig.pos), $(edge.dest.pos)), AT $inter\n")
                return edge, inter
            end
        end
        next_edge = nextedge(f, dir, edge)
		if counter >= 2 && !isnothing(commonvertex(edge, start_edge))
	        return nothing, [NaN, NaN]
        end
        edge = next_edge
		counter += 1
    end
    return
end


function oppositeface(f::Face, e::Edge; io=stdout)::Face
    #writenothing(io, "CHECKING OPPOSITE FACE OF: $(f)\n")
    #writenothing(io, "WITH MIRROR EDGE: $(e.id) ($(e.orig.pos), $(e.dest.pos))\n")
    #writenothing(io, "MIRROR RIGHT FACE: $(e.fr)\n")
    #writenothing(io, "MIRROR LEFT FACE: $(e.fl)\n")
    if e.fr == f
        #writenothing(io, "FOUND OPPOSITE: $(e.fl)\n")
        return e.fl
    elseif e.fl == f
        #writenothing(io, "FOUND OPPOSITE: $(e.fr)\n")
        return e.fr
    else
        throw("Face $(f) has a topology problem")
    end
    return
end


function joindcel(a::DCEL, b::DCEL)::DCEL
    append!(a.vertexlist, b.vertexlist)
    append!(a.edgelist, b.edgelist)
    append!(a.facelist, b.facelist)
    cleardcel!(b)
    return a
end



#---
function findedge(v::Vertex, f::Face, side::Symbol)::Edge
	edge = v.edge
    while edge.fr != f && edge.fl != f
		edge = nextpivot(edge, v, side)
    end
	return edge
end


function startedges(handler::Handler)
    infinite = isnothing(handler.current_vertex)
    if infinite
        return midpoint(handler.left_face, handler.right_face), handler.left_face.edge, handler.right_face.edge
    else
        return handler.current_vertex.pos, findedge(handler.left_vertex, handler.left_face, :ccw), findedge(handler.right_vertex, handler.right_face, :cw)
    end
    return
end


function highestintersection(D::DCEL, handler::Handler; io=stdout)::Tuple
    start, left_starter_edge, right_starter_edge = startedges(handler)
    perp_vector = perpvector(handler.right_face, handler.left_face)

	#writenothing(io, "\nCHECKING INTERSECTION OF DIAGRAMS WITH RAY\nSTART: $start\nVECTOR: $perp_vector\nIGNORE: $(handler.ignore)\n")
    #writenothing(io, "\nCHECKING LEFT DIAGRAM...\n")
    el, il = facerayintersection(handler.left_face, start, perp_vector, left_starter_edge, dir=:ccw, infinite=infinite, ignore=handler.ignore, io=io)
    #writenothing(io, "\nCHECKING RIGHT DIAGRAM...\n")
    er, ir = facerayintersection(handler.right_face, start, perp_vector, right_starter_edge, dir=:cw, infinite=infinite, ignore=handler.ignore, io=io)
    if isnothing(el) || ir[2] > il[2]
        split = createvertex!(D, ir)
		#writenothing(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT v$(split.id) $(ir)++\n")
        return true, split, splitedge!(D, er, split)
    elseif isnothing(er) || il[2] > ir[2]
        split = createvertex!(D, il)
        s = splitedge!(D, el, split)
		#writenothing(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT v$(split.id) $(il)++\n")
        return false, split, s
    else #sanity
        throw("")
    end
end


function killface!(face::Face, edge::Edge, dir::Symbol, stop_vertex::Vertex; io=stdout)::Nothing
    #writenothing(io, "KILLING FACE: $face\n")
    while true
        next_edge = nextedge(face, dir, edge)
        edge.dead = true
		#writenothing(io, "KILLED EDGE: e$(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        v = commonvertex(edge, next_edge)
        if stop_vertex in endpoints(edge)
            #writenothing(io, "STOPPED AT: $stop_vertex\n")
            return
        end
        v.dead = true
		#writenothing(io, "KILLED VERTEX: $v\n")
        edge = next_edge
    end
	#writenothing(io, "STOPPED AT: $stop_vertex\n")
    return
end


function updatehandler!(handler::Handler, edge::Edge, new_vertex::Vertex, 
    side::Bool, ignore::Tuple; io=stdout)::Nothing
    if side
        handler.right_face = oppositeface(handler.right_face, edge, io=io)
        handler.right_vertex = new_vertex
    else
        handler.left_face = oppositeface(handler.left_face, edge, io=io)
        handler.left_vertex = new_vertex
    end
    handler.current_vertex = new_vertex
    handler.current_joint = new_vertex.edge
    handler.ignore = ignore
    return
end

function openface(f::Face, dir::Symbol; io=stdout)::Tuple
    e = findframe(f)
    #writenothing(io, "OPENED: $e\n")
    next_edge = nextedge(f, dir, e)
    v = commonvertex(e, next_edge)
    if v == e.orig
        return e, :dest
    else
        return e, :orig
    end
    return
end


function updateray!(D::DCEL, ray::Edge, u::Vertex, left::Face, right::Face)::Nothing
	v = createdummyvertex!(D, u, perpangle(left, right))
	replacevertex!(ray.dest, v)
	settopology!(ray, orig=u, fl=left, fr=right)
	squeezeedge!(u, ray, update_faces=false)
	squeezeedge!(v, ray, update_faces=false)
	u.edge = v.edge = ray
	return
end


function weldboundary(t::Tuple, ray::Edge, side::Bool; io=stdout)::Nothing
	#writenothing(io, "WELDING $(t[1]) to $ray\n")
	v = ray.dest
	getfield(t[1], t[2]).dead = true
	unstickedge!(t[1], t[2])
	setfield!(t[1], t[2], v)
	if ray.ccwd == ray
		squeezeedge!(v, t[1], update_faces=false, previous=ray, next=ray)
	else
		if side
			squeezeedge!(v, t[1], update_faces=false, previous=ray, next=ray.ccwd)
		else
			squeezeedge!(v, t[1], update_faces=false, previous=ray.ccwd, next=ray)
		end
	end
	return
end


function switchvertex!(e::Edge, old::Vertex, new::Vertex)::Nothing
	if old == e.orig
		e.orig = new
	elseif old == e.dest
		e.dest = new
	else
		throw("")
	end
	return
end


function replacevertex!(old::Vertex, new::Vertex)::Nothing
	edges = vertexedges(old)
	for e in edges
		switchvertex!(e, old, new)
	end
	return
end


function cleardcel!(D::DCEL)::Nothing
	D.vertexlist = []
	D.edgelist = []
	D.facelist = []
	return
end



##---
function mergevoronoi(left::DCEL, right::DCEL, io)
    hl, hr, ll, lr = findextrema(left, right)
    #writenothing(io, "TOP FACES: ($hl, $hr)\n")
    #writenothing(io, "BOT FACES: ($ll, $lr)\n")

    D = joindcel(left, right)

    tl, tr, bl, br = openface.([hl, hr, ll, lr], [:ccw, :cw, :cw, :ccw], io=io)
    starter_left = getfield(tl[1], tl[2])
    starter_right = getfield(tr[1], tr[2])
    top_ray = createfloatingedge!(D)
    bot_ray = createfloatingedge!(D) #comment out if you cannot print the edgelist

    handler = Handler(hl, hr, starter_left, starter_right)
    right_first, split, s = highestintersection(D, handler, io=io)
    if right_first
        edge = cwface(hr, s[1]) == s[2] ? s[2] : s[1]
        killface!(handler.right_face, edge, :cw, starter_right, io=io)
        weldboundary(tr, top_ray, true, io=io)
    else
        edge = ccwface(hl, s[1]) == s[2] ? s[2] : s[1]
        killface!(handler.left_face, edge, :ccw, starter_left, io=io)
        weldboundary(tl, top_ray, false, io=io)
    end
    updateray!(D, top_ray, split, hl, hr)
    updatehandler!(handler, edge, split, right_first, s, io=io)
    handler.current_joint = top_ray
    # return D
    while handler.left_face != ll || handler.right_face != lr
        right_first, split, s = highestintersection(D, handler, io=io)
        joint = joinvertices!(D, handler.current_vertex, split, split_face=false, update_edges=false)
        split.edge = joint
        settopology!(joint, fr=handler.left_face, fl=handler.right_face)
        #writenothing(io, "CREATED VERTEX v$(split.id)\n")
        #writenothing(io, "CREATED JOINT $(joint.id): ($(joint.orig.pos), $(joint.dest.pos))\n")
        if right_first
            edge = joint.ccwd
            handler.right_face.edge = joint
            killface!(handler.right_face, edge, :cw, handler.right_vertex, io=io)
        else
            edge = joint.cwd
            handler.left_face.edge = joint
            killface!(handler.left_face, edge, :ccw, handler.left_vertex, io=io)
        end
        if handler.current_vertex == handler.right_vertex
            settopology!(joint, cwo=handler.current_joint)
            ccwset!(handler.current_joint, handler.current_vertex, joint)
        else
            settopology!(joint, ccwo=handler.current_joint)
            cwset!(handler.current_joint, handler.current_vertex, joint)
        end
        
        if right_first && handler.right_face == hr #This also works for the degenerate case!
            weldboundary(tr, top_ray, true, io=io)
        elseif !right_first && handler.left_face == hl 
            weldboundary(tl, top_ray, false, io=io)
        end
        
        updatehandler!(handler, edge, split, right_first, s, io=io)
    end
    weldboundary(bl, bot_ray, true, io=io)
    updateray!(D, bot_ray, handler.current_vertex, lr, ll)
    weldboundary(br, bot_ray, false, io=io)
    if handler.current_vertex == handler.right_vertex
        squeezeedge!(handler.current_vertex, bot_ray, update_faces=false, next=handler.current_joint.cwd, previous=handler.current_joint)
    else
        squeezeedge!(handler.current_vertex, bot_ray, update_faces=false, previous=handler.current_joint.ccwd, next=handler.current_joint)
    end

    hl.edge = tl[1]
    hr.edge = tr[1]
    ll.edge = bl[1]
    lr.edge = br[1]

    return D
end


function deleteghostvertex(v::Vertex)
    e1 = v.edge
    e2 = ccw(e1, v)

    if v == e2.orig
        u = e2.dest
        a = e2.ccwd
        b = e2.cwd
    else
        u = e2.orig
        a = e2.ccwo
        b = e2.cwo
    end

    unstickedge!(e2, :orig)
    unstickedge!(e2, :dest)

    v.dead = true
    e2.dead = true

    if v == e1.orig
        unstickedge!(e1, :orig)
        e1.orig = u
    else
        unstickedge!(e1, :dest)
        e1.dest = u
    end
    squeezeedge!(u, e1, previous=b, next=a, update_faces=false)
    if !e1.orig.original
        flipedge!(e1)
    end
    return
end



function deleteghostvertices(D::DCEL)
    t = findall(x -> getfield(x, :ghost), D.vertexlist)
    for i in t
        deleteghostvertex(D.vertexlist[i])
    end
    deleteat!(D.vertexlist, t)
    return
end


function flipedge!(e::Edge)::Nothing
    e.orig, e.dest = e.dest, e.orig
    e.cwo, e.cwd = e.cwd, e.cwo
    e.ccwo, e.ccwd = e.ccwd, e.ccwo
    e.fl, e.fr = e.fr, e.fl
    return
end


function voronoi(points::Array, io)
    if length(points) == 2
        return voronoitwopoints(points)
    elseif length(points) == 3
        return voronoithreepoints(points)
    else
        split = Int(floor(length(points)/2))
        points_left = points[1:split]
        points_right = points[split+1:end]
        vor_left = voronoi(points_left, io)
        vor_right = voronoi(points_right, io)
        #writenothing(io, "\n==MERGING $points_left AND $points_right==\n")
        vor = mergevoronoi(vor_left, vor_right, io)
        #=TODO: Is this faster than filtering after every deletion? 
        Or filtering right before returning the final diagram? I think it is=#
        filter!(x->!getfield(x, :dead), vor.edgelist) 
        filter!(x->!getfield(x, :dead), vor.vertexlist)
        return vor
    end
    return
end


function voronoihelper(points::Array; io=stdout)
    global TOTAL_EDGES, TOTAL_FACES, TOTAL_VERTICES
    TOTAL_EDGES = TOTAL_FACES = TOTAL_VERTICES = 0
    points = sort(points, by=x->x[1])
    out = voronoi(points, io)
    deleteghostvertices(out)
    filter!(x->!getfield(x, :dead), out.vertexlist)
    filter!(x->!getfield(x, :dead), out.edgelist)
    
    fixids!(out)
    return out
end


end
