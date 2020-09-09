module DStruct
import Base.show
using Plots, Statistics
export DCEL, Vertex, Edge, Face, Handler,

createvertex!, createdummyvertex!, splitedge!, addray!, joinvertices!, createfloatingedge!,#Constructors

squeezeedge!, settopology!, cwset!, ccwset!, #Topology

ccw, cw, ccwface, cwface, nextedge, #Face and edge traversal

faceedges, faceverticescw, faceverticesccw, vertexedges,#

mergeinfinitefaces!, joindcel, #

fixids!,

checkedge, checkdcel,

plotdcel,

facerayintersection,

midpoint, bisectorangles, distance,

findextrema, findstruts, findsupport, leftofline, findframe,

findside, getframesite,

perpangle, isframe, commonvertex, oppositeface, endpoints, unstickedge!,

voronoitwopoints, voronoithreepoints, #Elementary diagrams

openface, foobar, highestintersection, updatehandler!, weldboundary, updateray!, writenothing, cleardcel!,

averagepositions, sorted_ccw, uncommonvertices, sorted_ccw_1, sorted_ccw_2

#---
#Struct
mutable struct Vertex
    id::String
    pos::Array
    edge
    original::Bool
    dead::Bool
end
Vertex(id)::Vertex = Vertex(id, [NaN, NaN], nothing, true, false)
Vertex(id, pos)::Vertex = Vertex(id, pos, nothing, true, false)
Vertex(id, pos, edge)::Vertex = Vertex(id, pos, edge, true, false)
Vertex(id, pos, edge, original)::Vertex = Vertex(id, pos, edge, original, false)

mutable struct Face
    id::String
    edge
    site::Union{Array,Nothing}
end
Face(id)::Face = Face(id, nothing, nothing)
Face(id, edge)::Face = Face(id, edge, nothing)


mutable struct Edge
    id::String
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
    edgelist::Array
    vertexlist::Array
    facelist::Array
end
DCEL() = DCEL([], [], [])


mutable struct Handler
    left_face::Union{Face,Nothing}
    right_face::Union{Face,Nothing}
    left_vertex::Union{Vertex,Nothing}
    right_vertex::Union{Vertex,Nothing}
    current_vertex::Union{Vertex,Nothing}
    current_joint::Union{Edge,Nothing}
    side::Union{Symbol,Nothing}
    ignore::Array
end
Handler(lf::Face, rf::Face, lv::Vertex, rv::Vertex)::Handler = Handler(lf, rf, lv, rv, nothing, nothing, nothing, [])



#---
#Printing
function printnothing(e::Any)::String
    p = String
    if e isa Edge
        p = "e"
    elseif e isa Face
        p = "f"
    elseif e isa Vertex
        p = "v"
    end
    return isnothing(e) ? "nothing" : "$(p)$(e.id)"
end


function show(io::IO, v::Vertex)::Nothing
	    print(io, "v$(v.id): e=", printnothing(v.edge))
    print(io, " pos=$(v.pos) orig=$(v.original)")
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
	print(io, " edge=e$(f.edge.id)")
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
@inline function cross2d(u::T, v::T) where T<:Union{Tuple, Array}
    return u[1]*v[2] - u[2]*v[1]
end


@inline function dot(u::Array, v::Array)::Float64
    return u[1]*v[1] + u[2]*v[2]
end


@inline function norm(a::Array)::Float64
    return hypot(a...)
end


@inline function distance(p::Array, q::Array)::Float64
	return norm(p .-q)
end


@inline function distance(p::Vertex, q::Vertex)::Float64
	return distance(p.pos, q.pos)
end


@inline function cosatan(y::Number, x::Number)::Float64
    return x/hypot(x,y)
end


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


function pointccw(array::Array)::Bool
    sum = 0
    for i in 2:length(array)
        sum += (array[i][1] - array[i-1][1])*(array[i][2] + array[i-1][2])
    end
    sum += (array[1][1] - array[end][1])*(array[1][2] + array[end][2])
    return sum <= 0
end


function midpoints(u::Array, v::Array, w::Array)::Tuple
    return mean([u,v]), mean([v,w]), mean([w,u])
end


function angleccw(a::Number, b::Number, c::Number)::Bool
    return sin(a-b) + sin(b-c) + sin(c-a)  <= 0
end


function circletwopointsradius(p::Array, q::Array, r::Number)
    x1, y1 = p
    x2, y2 = q
	q = hypot(x2-x1, y2-y1)
    y3 = (y1+y2)/2
    x3 = (x1+x2)/2
    basex = sqrt(r^2 - (q/2)^2) * (y1-y2)/q
    basey = sqrt(r^2 - (q/2)^2) * (x2-x1)/q
    return [x3+basex, y3+basey], [x3-basex, y3-basey]
end


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
function createvertex!(D::DCEL, pos::Array, original::Bool=true)::Vertex
     new_vertex = Vertex("$(length(D.vertexlist)+1)", pos, nothing, original)
    push!(D.vertexlist, new_vertex)
    return new_vertex
end


function createedge!(D::DCEL, u::Vertex, v::Vertex)::Edge
    new_edge = Edge("$(length(D.edgelist)+1)", u, v)
    settopology!(new_edge)
    push!(D.edgelist, new_edge)
    return new_edge
end


function createfloatingedge!(D::DCEL)::Edge
	u = Vertex("x")
	v = Vertex("y")
	e = joinvertices!(D, u, v)
	settopology!(e, cwd=e, ccwd=e, cwo=e, ccwo=e)
	return e
end


function createface!(D::DCEL)::Face
    new_face = Face("$(length(D.facelist)+1)")
    push!(D.facelist, new_face)
    return new_face
end


function createdummyvertex!(D::DCEL, u::Vertex, angle::Number)::Vertex
    vector = u.pos .+ [cos(angle),sin(angle)]
    v = createvertex!(D, vector, false)
    return v
end


function fixids!(D::DCEL)::Nothing
    for i in 1:length(D.vertexlist)
        D.vertexlist[i].id = "$(i)"
    end
    for i in 1:length(D.edgelist)
        D.edgelist[i].id = "$(i)"
    end
    for i in 1:length(D.facelist)
        D.facelist[i].id = "$(i)"
    end
    return
end


function averagepositions(D::DCEL)
	pos = getfield.(D.vertexlist, :pos)
	mean_pos = mean(pos)
	dist = [abs.(p - mean_pos) for p in pos]
	mean_dist = mean(dist)
	return mean_pos, mean_dist
end


function plotdcel(D::DCEL; faces::Bool=false, dead_edges::Bool=false,
				  dead_vertices::Bool=false, bounds::Bool=true, sites::Bool=false,
				  normalize::Bool=true, scale::Number=1, font_size::Int=12,
				  labels::Bool=true, line::Symbol=:arrow, kwargs...)::Plots.Plot

	avg_pos, avg_dist = averagepositions(D)
	xlims = avg_pos[1] .+ [-1, 1] .*avg_dist
	ylims = avg_pos[2] .+ [-1, 1] .*avg_dist
	p = plot(leg=false, xlims=xlims, ylims=ylims; kwargs...)


	edge_pos_x = Array{Float64}(undef,0)
	edge_pos_y = Array{Float64}(undef,0)
	edge_color = []
	edge_text = []
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
		text = labels ? (ave[1], ave[2], "e$(e.id)", Plots.font("Sans",font_size)) : []
        if (dead_edges || !e.dead) && (bounds || !bound)
			push!(edge_pos_x, start[1], finish[1], NaN)
			push!(edge_pos_y, start[2], finish[2], NaN)
			push!(edge_color, color)
			push!(edge_text, text)
        end
    end
	plot!(p, edge_pos_x, edge_pos_y,color=edge_color, line=line)

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
		text = labels ? [Plots.text("\nv$(v.id)",font_size)] : []
		if (dead_vertices || !v.dead) && (bounds || v.original)
			push!(vertex_pos_x, pos[1])
			push!(vertex_pos_y, pos[2])
			push!(vertex_col, color)
			push!(vertex_text, text)
		end
	end
	scatter!(p, vertex_pos_x, vertex_pos_y, color=vertex_col)

	face_pos_x = []
	face_pos_y = []
	face_text = []
	site_pos_x = []
	site_pos_y = []
    for f in D.facelist
        if faces
            vertices = faceverticesccw(f)
            x, y = mean(getfield.(vertices, :pos))
			push!(face_pos_x, x)
			push!(face_pos_y, y)
			push!(face_text, "f$(f.id)")
        end
        if sites && !isnothing(f.site)
			push!(site_pos_x, f.site[1])
			push!(site_pos_y, f.site[2])
        end
    end
	annotate!(p, face_pos_x, face_pos_y, face_text)
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
		throw("")
	end
end


function ccw(e::Edge, v::Vertex)::Edge
	if v == e.orig
		return e.ccwo
	elseif v == e.dest
		return e.ccwd
	else #sanity
		throw("")
	end
end


function nextpivot(e::Edge, v::Vertex, side::Symbol)::Edge
	@assert side in [:cw,:ccw]
	return getfield(Main, side)(e, v)
end


function cwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? settopology!(e, cwo=new_edge) : settopology!(e, cwd=new_edge)
    return
end


function ccwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? settopology!(e, ccwo=new_edge) : settopology!(e, ccwd=new_edge)
    return
end


@inline function endpoints(u::Edge)::Array
    return [u.orig, u.dest]
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


function sorted_ccw_2(u::Edge, v::Edge, w::Edge)::Bool
    a, b, c, p = getfield.(uncommonvertices(u, v, w), :pos)
    a = a - p
    b = b - p
    c = c - p
    return cross2d(a, b) * sqrt(c[1]^2+c[2]^2) + cross2d(b, c) * sqrt(a[1]^2+a[2]^2) + cross2d(c, a) * sqrt(b[1]^2+b[2]^2) >= 0
end


function sorted_ccw(u::Edge, v::Edge, w::Edge)::Bool #TODO: Compare
    return sorted_ccw_2(u, v, w)
end

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


function getstrut(v::Vertex)::Edge #TODO: add sanity
    @assert !v.original "v is not a dummy vertex"
    edge = v.edge
    while !isstrut(edge)
        edge = ccw(edge, v)
    end
    return edge
end


function findframe(f::Face; dir::Symbol=:ccw)::Edge
	edge = f.edge
	while !isframe(edge)
		edge = nextedge(f, dir, edge)
	end
	return edge
end


function findstruts(f::Face; dir::Symbol=:ccw)::Tuple #TODO: add sanity
    edge = findframe(f)
    return cwface(f, edge), ccwface(f, edge)
end


function normalizedummy(v::Vertex)::Array
    @assert !v.original "v is not a dummy vertex"
    e = getstrut(v)
    u, v = e.orig.original ? (e.orig, e.dest) : (e.dest, e.orig)
    return v.pos - u.pos
end


vertexccw(vertices::Array) = pointccw(getfield.(vertices, :pos))


function isframe(e::Edge)::Bool
    # return all(.!endfield(e, :original))
	return !e.orig.original  && !e.dest.original
end


function isstrut(e::Edge)::Bool
    # a = endfield(e, :original)
    # return xor(a[1], a[2])
	return (e.orig.original && !e.dest.original) || (!e.orig.original && e.dest.original)
end


#Magic
function squeezeedge!(v::Vertex, e::Edge, update_faces::Bool=true;
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


function add_join_vertex!(D::DCEL, new::Array, old_vertex::Union{Vertex, Nothing}=nothing)::Vertex
    new_vertex = createvertex!(D, new)

    if isnothing(old_vertex)
        return new_vertex
    else
        @assert old_vertex in D.vertexlist "Given vertex is not part of DCEL."
    end

    new_edge = createedge!(D, old_vertex, new_vertex)
    new_vertex.edge = new_edge

    squeezeedge!(old_vertex, new_edge)
    if isnothing(old_vertex.edge)
        f = Face("i", new_edge)
        push!(D.facelist, f)
        new_edge.fr = new_edge.fl = f
    end
    old_vertex.edge = new_edge
    return new_vertex
end


function endfield(e::Edge, f::Symbol)
    return getfield.(endpoints(e), f)
end


edgelength(e::Edge) = norm(mean(endfield(e, :pos)))::Float64


edgemidpoint(e::Edge) = mean(endfield(e, :pos))::Array


function joinvertices!(D::DCEL, u::Vertex, v::Vertex; po::Union{Edge,Nothing}=nothing, pd::Union{Edge,Nothing}=nothing, no::Union{Edge,Nothing}=nothing, nd::Union{Edge,Nothing}=nothing, split_face::Bool=true, update_edges::Bool=true)
    @assert u != v "Cannot join a vertex to itself"
    new_edge = createedge!(D, u, v)
    squeezeedge!(u, new_edge, previous=po, next=no)
    squeezeedge!(v, new_edge, previous=pd, next=nd)
    disconnected = isnothing(u.edge) || isnothing(v.edge)
	if update_edges
    	u.edge = v.edge = new_edge
	end
    if disconnected
        return new_edge
    end

    if split_face
        f = new_edge.fr
        f1 = Face("$(length(D.facelist)+1)", new_edge)
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

        f2 = Face("$(length(D.facelist)+2)", new_edge)
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

        filter!(x->x != f, D.facelist)
        push!(D.facelist, f1, f2)
    end
    return new_edge
end


function cwface(f::Face, edge::Union{Edge,Nothing}=nothing)::Edge
    if isnothing(edge)
        return f.edge
    end
	# println("CURRENTLY AT $edge")
    if edge.fr == edge.fl
        return edge.ccwo != edge ? edge.ccwo : edge.ccwd
    elseif f == edge.fr
		# println("JUMPED TO: $(edge.ccwd)")
        return edge.ccwd
    elseif f == edge.fl
		# println("JUMPED TO: $(edge.ccwo)")
        return edge.ccwo
    else #sanity
        throw("$f has a topology problem with edge $edge")
    end
    return
end


function ccwface(f::Face, edge::Union{Edge,Nothing}=nothing)::Edge
    if isnothing(edge)
        return f.edge
    end
    if edge.fr == edge.fl
        return edge.cwo != edge ? edge.cwo : edge.cwd
    elseif f == edge.fr
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

function firstlastdummy(D::DCEL)::Tuple
    v1 = findfirst(x -> !getfield(x, :original), D.vertexlist)
    v2 = findlast(x -> !getfield(x, :original), D.vertexlist)
    return isnothing(v1) ? (nothing, nothing) : (D.vertexlist[v1], D.vertexlist[v2])
end

#If site data is not available
# function findinfiniteface(D::DCEL)::Union{Face,Nothing}
#     v1, v2 = firstlastdummy(D)
#     if isnothing(v1)
#         return nothing
#     end
#     if v1 == v2
#         return v1.edge.fr
#     else
#         v = v1
#         edge = v.edge
#         while !isframe(edge) #TODO: add loop protection
#             edge = ccw(edge, v)
#         end
#         fr = edge.fr
#         fl = edge.fl
#         if isframe(ccwface(fr, edge))
#             return fr
#         elseif isframe(ccwface(fl, edge))
#             return fl
#         end
#     end
#     return
# end


function findinfiniteface(D::DCEL)::Face
	index = findfirst(x->isnothing(getfield(x, :site)), D.facelist)
	return D.facelist[index]
end


function deleteedge!(D::DCEL, e::Edge, merge_faces::Bool=true)::Nothing
    if merge_faces
        f = Face("$(length(D.facelist))")
        f.edge = e.ccwo
        resetfacelist!(D, f, e.fl)
        resetfacelist!(D, f, e.fr)
        push!(D.facelist, f)
    end
    ccwset!(e.cwo, e.orig, e.ccwo)
    cwset!(e.ccwo, e.orig, e.cwo)
    ccwset!(e.cwd, e.dest, e.ccwd)
    cwset!(e.ccwd, e.dest, e.cwd)
    filter!(x->x != e, D.edgelist)
    return
end


function resetface(face::Face, new_face::Face)::Nothing
    for e in faceedges(face, :ccw)
        if e.fr == face
            e.fr = new_face
        else
            e.fl = new_face
        end
    end
    return
end


function resetfacelist!(D::DCEL, face::Face, new_face::Face)::Nothing
    resetface(face, new_face)
    filter!(x->x !=face, D.facelist)
    return
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


function addray!(D::DCEL, u::Vertex, angle::Number, f::Union{Face,Nothing}=nothing, split_face::Bool=true)::Edge
    v1, v2 = firstlastdummy(D)
    v = createdummyvertex!(D, u, angle)
    new_edge = nothing
    if v1 == v2
        new_edge = joinvertices!(D, u, v)
        if !isnothing(v1)
            outer = joinvertices!(D, v1, v)
            joinvertices!(D, v, v1, po=outer, pd=getstrut(v1))
        end
    else
        outeredge = nothing
        if isnothing(f)
            edge = u.edge
            outeredge = checkedgefaces(edge)
            while isnothing(outeredge)
                edge = ccw(edge, u)
                outeredge = checkedgefaces(edge)
            end
            f = isframe(ccwface(outeredge.fr, outeredge)) ? outeredge.fl : outeredge.fr
        else
            outeredge = findfaceframe(f)
        end
        e1, e2 = splitedge!(D, outeredge, v)
        previous = ccwface(f, e1) == e2 ? e2 : e1
        new_edge = joinvertices!(D, u, v, split_face=split_face, pd=previous)
    end
    return new_edge
end


function vertexinface(v::Vertex, f::Face)::Bool
    edge = f.edge
    while !(v in endpoints(edge))
        edge = ccwface(f, edge)
        if edge == f.edge
            return false
        end
    end
    return true
end


function checkedgefaces(e::Edge)::Union{Face,Nothing}
    a = findfaceframe(e.fl, e)
    b = findfaceframe(e.fr, e)
    if !isnothing(a)
        return a
    elseif !isnothing(b)
        return b
    else
        return nothing
    end
end


function splitedge!(D::DCEL, e::Edge, split::Vertex)::Array
    new_edge = Edge("$(e.id)b", split, e.dest)
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

    e.id = "$(e.id)a"
    e.dest.edge = new_edge
    e.dest = split
    e.cwd = e.ccwd = new_edge
    return [e, new_edge]
end


function line(a::Array, b::Array, c::Array)::Bool
    return ((b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1])) > 0
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
	if v.edge == e
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


function voronoitwopoints(points::Array, p::Int=2)::DCEL
    u, v = points
    D = DCEL()
    m = createvertex!(D, mean([u,v]))
    angle = atan(v[2]-u[2], v[1]-u[1]) + pi/2
    r = addray!(D, m, angle)
    addray!(D, m, angle+pi)
	r.fl.site = u
	r.fr.site = v
    return D
end


function voronoithreepoints(points::Array)::DCEL
    u, v, w = points
    a, b, c = pointccw(points) ? bisectorangles(v, u, w) : bisectorangles(u, v, w) #somehow this works
    D = DCEL()
    m = createvertex!(D, circlethreepoints(u, v, w))
    p = addray!(D, m, a)
    r = addray!(D, m, c)
    q = addray!(D, m, b, angleccw(a,b,c) ? p.fl : p.fr)
    if pointccw(points)
        p.fr.site = w
        q.fr.site = v
        r.fr.site = u
    else
        p.fl.site = v
        q.fl.site = w
        r.fl.site = u
    end
    return D
end


function getframeface(e::Edge)::Face
    if isnothing(e.fr.site)
        return e.fl
    elseif isnothing(e.fl.site)
        return e.fr
    else #sanity
        throw("")
    end
    return
end


function getframesite(e::Edge)::Array
    return getframeface(e).site
end


function findside(D::DCEL, max::Bool)
    frames = findall(x-> isframe(x), D.edgelist)
    edges = D.edgelist[frames]
    f = max ? argmax : argmin
	return edges[f(hcat(getfield.(getframeface.(edges), :site)...)[1,:])]
end


function findsupport(left::DCEL, right::DCEL, f::Face, top::Bool)::Tuple#TODO: Optimize
	edge_l = findside(left, true)
    edge_r = findside(right, false)
    f_l = top ? cwface : ccwface
    f_r = top ? ccwface : cwface
    next_l = f_l(f, edge_l)
    next_r = f_r(f, edge_r)
    cond = top ? !leftofline : leftofline
    while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_l)) || cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_r))
        while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_l))
            edge_l = next_l
            next_l = f_l(f, edge_l)
        end
        while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_r))
            edge_r = next_r
            next_r = f_r(f, edge_r)
        end
    end
    return getframeface(edge_l), getframeface(edge_r)
end


function findextrema(left::DCEL, right::DCEL, f::Face)::Tuple
    return (findsupport(left, right, f, true)..., findsupport(left, right, f, false)...)
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


# function edgerayintersection(edge::Edge, q::Array, angle::Float64, infinite::Bool=false; io=stdout)::Array
#     p = edge.orig.pos
#     r = edge.dest.pos - p
#     s = [cos(angle), sin(angle)]
#     num = q - p
#     den = cross2d(r, s)
#     t = cross2d(num, s)/den
#     u = cross2d(num, r)/den
# 	cond = isstrut(edge) ? 0<=t : 0<=t<=1
# 	return (u>=0 || infinite) && cond ? p+t*r : [NaN,NaN]
# end


# function facerayintersection(f::Face, start::Array, angle::Number, edge::Edge;
# 	dir::Symbol=:ccw, infinite::Bool=false, ignore::Array=[], io=stdout)::Tuple
# 	start_vertices = endpoints(edge)
# 	start_edge = edge
# 	counter = 0
#     #writenothing(io, "CHECKING INTERSECTION OF FACE $(f) WITH IGNORE $(getfield.(ignore, :id))\n")
#     while true
#         #writenothing(io, "CHECKING EDGE $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
#         if !edge.dead && !isframe(edge) && !(edge in ignore)
#             inter = edgerayintersection(edge, start, angle, infinite, io=io)
# 			#writenothing(io, "INTERSECTION: $inter\n")
#             if !isnan(inter[1])
#                 #writenothing(io, "INTERSECTED $(edge.id) ($(edge.orig.pos), $(edge.dest.pos)), AT $inter\n")
#                 return edge, inter
#             end
#         end
#         next_edge = nextedge(f, dir, edge)
# 		if counter >= 2 && !isnothing(commonvertex(edge, start_edge))
# 	        return nothing, [NaN, NaN]
# 		end

#         edge = next_edge
# 		counter += 1
#     end
#     return
# end


function edgerayintersection(edge::Edge, q::Array, vector::Array, infinite::Bool=false; io=stdout)::Array
    p = edge.orig.pos
    r = edge.dest.pos - p
    num = q - p
    den = cross2d(r, vector)
    t = cross2d(num, vector)/den
    u = cross2d(num, r)/den
	cond = isstrut(edge) ? 0<=t : 0<=t<=1
	return (u>=0 || infinite) && cond ? p+t*r : [NaN,NaN]
end


function facerayintersection(f::Face, start::Array, vector::Array, edge::Edge;
	dir::Symbol=:ccw, infinite::Bool=false, ignore::Array=[], io=stdout)::Tuple
	start_vertices = endpoints(edge)
	start_edge = edge
	counter = 0
    #writenothing(io, "CHECKING INTERSECTION OF FACE $(f) WITH IGNORE $(getfield.(ignore, :id))\n")
    while true
        #writenothing(io, "CHECKING EDGE $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
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
    if e.fr.site == f.site
        #writenothing(io, "FOUND OPPOSITE: $(e.fl)\n")
        return e.fl
    elseif e.fl.site == f.site
        #writenothing(io, "FOUND OPPOSITE: $(e.fr)\n")
        return e.fr
    else
        throw("Face $(f) has a topology problem")
    end
    return
end


function mergeinfinitefaces!(a::DCEL, b::DCEL)::Face
    fi = Face("i")
    fa = findinfiniteface(a)
    fb = findinfiniteface(b)
    resetfacelist!(a, fa, fi)
    resetfacelist!(b, fb, fi)
    fi.edge = fa.edge
    return fi
end


function joindcel(a::DCEL, b::DCEL)::DCEL
    vertexlist = vcat(a.vertexlist, b.vertexlist)
    edgelist = vcat(a.edgelist, b.edgelist)
    facelist = vcat(a.facelist, b.facelist)
    return DCEL(edgelist, vertexlist, facelist)
end



#---
function findedge(v::Vertex, f::Face, side::Symbol)::Edge
	edge = v.edge
	while edge.fr != f && edge.fl != f
		edge = nextpivot(edge, v, side)
	end
	return edge
end


function highestintersection(D::DCEL, handler::Handler; io=stdout)::Tuple
	infinite = isnothing(handler.current_vertex)
	if infinite
        start = midpoint(handler.left_face, handler.right_face)
		left_starter_edge = handler.left_face.edge
		right_starter_edge = handler.right_face.edge
    else
        start = handler.current_vertex.pos
		left_starter_edge = findedge(handler.left_vertex, handler.left_face, :ccw)
		right_starter_edge = findedge(handler.right_vertex, handler.right_face, :cw)
    end
    # angle = perpangle(handler.right_face, handler.left_face)
    perp_vector = perpvector(handler.right_face, handler.left_face)
	#writenothing(io, "\nCHECKING INTERSECTION OF DIAGRAMS WITH RAY\nSTART: $start\nANGLE: $(rad2deg(angle))\nIGNORE: $(handler.ignore)\n")
    #writenothing(io, "\nCHECKING LEFT DIAGRAM...\n")
    el, il = facerayintersection(handler.left_face, start, perp_vector, left_starter_edge, dir=:ccw, infinite=infinite, ignore=handler.ignore, io=io)
    #writenothing(io, "\nCHECKING RIGHT DIAGRAM...\n")
    er, ir = facerayintersection(handler.right_face, start, perp_vector, right_starter_edge, dir=:cw, infinite=infinite, ignore=handler.ignore, io=io)
    dr = distance(start, ir)
    dl = distance(start, il)
    s = nothing
    if all(isnan.(il)) || ir[2]>il[2]
        starter_right = :right
        split = createvertex!(D, ir)
		#writenothing(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT v$(split.id) $(ir)++\n")
        s = splitedge!(D, er, split)
    elseif all(isnan.(ir)) || il[2]>ir[2]
        starter_right = :left
        split = createvertex!(D, il)
		#writenothing(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT v$(split.id) $(il)++\n")
        s = splitedge!(D, el, split)
    else #sanity
        throw("")
    end
    return starter_right, split, s
end


function killface!(face::Face, edge::Edge, dir::Symbol, stop_vertex::Union{Vertex}; io=stdout)::Nothing
    #writenothing(io, "KILLING FACE: $face\n")
    while true
        next_edge = dir==:cw ? cwface(face, edge) : ccwface(face, edge)
        edge.dead = true
		#writenothing(io, "KILLED EDGE: e$(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        v = commonvertex(edge, next_edge)
        if stop_vertex in endpoints(edge)
			#writenothing(io, "STOPPED AT: $stop_vertex\n")
            return
        end
		if isnothing(v)
			throw("$edge, $next_edge don't have a common vertex!")
		end
        v.dead = true
		#writenothing(io, "KILLED VERTEX: $v\n")
        edge = next_edge
    end
	#writenothing(io, "STOPPED AT: $stop_vertex\n")
    return
end


function updatehandler!(handler::Handler, edge::Edge, new_vertex::Vertex, side::Symbol,
	stop_left::Vertex, stop_right::Vertex; io=stdout)::Nothing
    if side == :right
        killface!(handler.right_face, edge, :cw, stop_right, io=io)
        handler.right_face = oppositeface(handler.right_face, edge, io=io)
        handler.right_vertex = new_vertex
    elseif side == :left
        killface!(handler.left_face, edge, :ccw, stop_left, io=io)
        handler.left_face = oppositeface(handler.left_face, edge, io=io)
        handler.left_vertex = new_vertex
    else
        throw("")
    end
    handler.current_vertex = new_vertex
    handler.side = side
    handler.current_joint = new_vertex.edge
    return
end


function foobar(D::DCEL, handler::Handler;
	stop_left::Union{Vertex,Nothing} = handler.left_vertex,
	stop_right::Union{Vertex,Nothing} = handler.right_vertex, io=stdout)::Nothing
	right_first, split, s = highestintersection(D, handler, io=io)
    joint = joinvertices!(D, handler.current_vertex, split, split_face=false, update_edges=false)
	split.edge = joint
    settopology!(joint, fr=handler.left_face, fl=handler.right_face)
	#writenothing(io, "CREATED VERTEX v$(split.id)\n")
    #writenothing(io, "CREATED JOINT $(joint.id): ($(joint.orig.pos), $(joint.dest.pos))\n")
    if right_first == :right
        handler.right_face.edge = joint
    elseif right_first == :left
        handler.left_face.edge = joint
    end
    edge = getfield(joint, right_first==:right ? :ccwd : :cwd)
	if handler.side == :right
        joint.cwo = handler.current_joint
        ccwset!(handler.current_joint, handler.current_vertex, joint)
    else
        joint.ccwo = handler.current_joint
        cwset!(handler.current_joint, handler.current_vertex, joint)
    end
    updatehandler!(handler, edge, split, right_first, stop_left, stop_right, io=io)
    handler.ignore = s
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
	v = createdummyvertex!(D, u, perpangle(left,right))
	replacevertex!(ray.dest, v)
	settopology!(ray, orig=u, fl=left, fr=right)
	squeezeedge!(u, ray, false)
	squeezeedge!(v, ray, false)
	u.edge = v.edge = ray
	return
end


function weldboundary(D::DCEL, t::Tuple, ray::Edge, side::Symbol; io=stdout)::Nothing
	#writenothing(io, "WELDING $(t[1]) to $ray\n")
	v = ray.dest
	getfield(t[1], t[2]).dead = true
	unstickedge!(t[1], t[2])
	setfield!(t[1], t[2], v)
	if ray.ccwd == ray
		squeezeedge!(v, t[1], false, previous=ray, next=ray)
	else
		if side == :right
			squeezeedge!(v, t[1], false, previous=ray, next=ray.ccwd)
		else
			squeezeedge!(v, t[1], false, previous=ray.ccwd, next=ray)
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


function writenothing(io, s::String)::Nothing
	if !isnothing(io)
		write(io, s)
	end
	return
end


function cleardcel!(D::DCEL)::Nothing
	D.vertexlist = []
	D.edgelist = []
	D.facelist = []
	return
end

end
