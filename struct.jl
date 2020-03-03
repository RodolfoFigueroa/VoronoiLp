module DStruct
import Base.show
using Plots, Statistics
export DCEL, Vertex, Edge, Face,

createvertex!, splitedge!, addray!, joinvertices!, #Constructors

squeezeedge!, settopology!, #Topology

ccw, cw, ccwface, cwface, #Face and edge traversal

faceedgesccw, faceedgescw, faceverticescw, faceverticesccw, #

mergeinfinitefaces!, joindcel, #

fixids!, checkdcel, plotdcel,

facerayintersection, rayrayintersection, infinitefaceintersection,

midpoint, bisectorangles,

findextrema, findstruts, findsupport, leftofline,

perpangle, hideedge!, isframe, commonvertex, oppositeface, endpoints, unstickedge!,

voronoitwopoints, voronoithreepoints

##
#Struct
mutable struct Vertex
    id::String
    pos::Array
    edge
    original::Bool
    dead::Bool
end
Vertex(id) = Vertex(id, [NaN, NaN], nothing, true, false)
Vertex(id, pos) = Vertex(id, pos, nothing, true, false)
Vertex(id, pos, edge) = Vertex(id, pos, edge, true, false)
Vertex(id, pos, edge, original) = Vertex(id, pos, edge, original, false)

mutable struct Face
    id::String
    edge
    site::Union{Array,Nothing}
end
Face(id) = Face(id, nothing, nothing)
Face(id, edge) = Face(id, edge, nothing)

mutable struct Edge
    id::String
    orig::Vertex
    dest::Vertex
    cwo#::Union{Edge,Nothing}
    ccwo#::Union{Edge,Nothing}
    cwd#::Union{Edge,Nothing}
    ccwd#::Union{Edge,Nothing}
    fr::Union{Face,Nothing}
    fl::Union{Face,Nothing}
    dead::Bool
end
Edge(id, orig, dest) = Edge(id, orig, dest, nothing, nothing, nothing, nothing, nothing, nothing, false)

mutable struct DCEL
    edgelist::Array
    vertexlist::Array
    facelist::Array
end
DCEL() = DCEL([], [], [])

##
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

##
#Geometry
function cross2d(u::Array, v::Array)::Float64
    return u[1]*v[2] - u[2]*v[1]
end

function dot(u::Array, v::Array)::Float64
    return u[1]*v[1] + u[2]*v[2]
end

function norm(a::Array)::Float64
    return hypot(a...)
end

function distance(p::Array, q::Array)::Float64
	return norm(p .-q)
end

function cosatan(y::Number, x::Number)::Float64
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
    return sin(a-b) - sin(a-c) + sin(b-c) <= 0
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

function circlethreepoints(p::Array, q::Array, r::Array)
    x1,y1 = p
    x2,y2 = q
    x3,y3 = r
    den = 2*(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))
    x = (x1^2+y1^2)*(y3-y2)+(x2^2+y2^2)*(y1-y3)+(x3^2+y3^2)*(y2-y1)
    y = (x1^2+y1^2)*(x2-x3)+(x2^2+y2^2)*(x3-x1)+(x3^2+y3^2)*(x1-x2)
    return [x/den,y/den]
end

function circlethreepoints(p::Face, q::Face, r::Face)
	return circlethreepoints(p.site, q.site, r.site)
end

##
#
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

function plotdcel(D::DCEL; faces::Bool=false, dead_edges::Bool=false, dead_vertices::Bool=false, bounds::Bool=true, ratio::Any=:equal, sites::Bool=false)
    p = plot(leg=false)
    for v in D.vertexlist
        color = v.dead ? :gray : v.original ? :black : :red
        if dead_vertices || !v.dead
            scatter!(p, [v.pos[1]], [v.pos[2]], series_annotations=[Plots.text("\nv$(v.id)"), :bottom], color=color, aspect_ratio=ratio)
        end
    end
    for e in D.edgelist
        bound = isframe(e)
        color = e.dead ? :gray : bound ? :red : :black
        ave = edgecenter(e)
        if (dead_edges || !e.dead) && (bounds || !bound)
            plot!(p, [e.orig.pos[1],e.dest.pos[1]], [e.orig.pos[2],e.dest.pos[2]], line=:arrow, annotations = (ave[1], ave[2], "e$(e.id)"), linecolor=color, aspect_ratio=ratio)
        end
    end
    for f in D.facelist
        if faces
            vertices = faceverticesccw(f)
            x,y = mean(getfield.(vertices, :pos))
            annotate!([x],[y], "f$(f.id)")
        end
        if sites && !isnothing(f.site)
            scatter!([f.site[1]], [f.site[2]], color=:blue)
        end
    end
    return p
end

function checkdcel(D::DCEL)
    for v in D.vertexlist
        @assert v in endpoints(v.edge)
    end
    for e in D.edgelist
        @assert cw(e.ccwo, e.orig) == e "ccwo->cwo $(e)"
        @assert ccw(e.cwo, e.orig) == e "cwo->ccwo $(e)"
        @assert cw(e.ccwd, e.dest) == e "ccwd->cwd $(e)"
        @assert ccw(e.cwd, e.dest) == e "cwd->ccwd $(e)"
        @assert e.orig != e.dest
        if isstrut(e)
            @assert !e.dest.original
			@assert e.orig.original
        end
    end
    for f in D.facelist
        edgesccw = faceedgesccw(f)
        edgescw = faceedgescw(f)
        @assert isempty(setdiff(edgescw, edgesccw))
        for e in edgescw
            @assert f==e.fr || f==e.fl "$(f):$(e)"
        end
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

function cwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? e.cwo=new_edge : e.cwd=new_edge
    return
end

function ccwset!(e::Edge, v::Vertex, new_edge::Edge)::Nothing
    v==e.orig ? e.ccwo=new_edge : e.ccwd=new_edge
    return
end

function endpoints(u::Edge)::Set
    return Set([u.orig, u.dest])
end

function edgeangle(u::Edge)::Float64
    return atan(u.dest.pos[2]-u.orig.pos[2], u.dest.pos[1]-u.orig.pos[1])
end

function pivotangle(u::Edge, p::Vertex)::Float64
    return edgeangle(u) + (p==u.orig ? 0 : pi)
end

function returnempty(s::Set)
    return isempty(s) ? nothing : pop!(s)
end

function commonvertex(u::Edge, edges...)::Union{Vertex,Nothing}
    f = endpoints(u)
    for e in endpoints.(edges)
        intersect!(f, e)
    end
    return returnempty(f)
end

function uncommonvertices(u::Edge, v::Edge, w::Edge)
    s1 = endpoints(u)
    s2 = endpoints(v)
    s3 = endpoints(w)
    s = intersect(s1, s2, s3)
    return returnempty(setdiff(s1,s)), returnempty(setdiff(s2,s)), returnempty(setdiff(s3,s)), pop!(s)
end

function sorted_ccw(u::Edge, v::Edge, w::Edge)::Bool
    pivot = commonvertex(u, v, w)
    a = pivotangle(u, pivot)
    b = pivotangle(v, pivot)
    c = pivotangle(w, pivot)
    return sin(a-b) - sin(a-c) + sin(b-c) <= 0
end

# function sorted_ccw2(u::Edge, v::Edge, w::Edge)::Bool
#     a, b, c, d = uncommonvertices(u, v, w)
#     return (b.pos[1]-a.pos[1])*(b.pos[2]+a.pos[2])+(c.pos[1]-b.pos[1])*(c.pos[2]+b.pos[2])+(d.pos[1]-c.pos[1])*(d.pos[2]+c.pos[2])+(a.pos[1]-d.pos[1])*(a.pos[2]+d.pos[2]) <= 0
# end
#
# function sorted_ccw3(u::Edge, v::Edge, w::Edge)::Bool
#     p,q,r,z = uncommonvertices(u, v, w)
#     A = p.pos-z.pos
#     B = q.pos-z.pos
#     C = r.pos-z.pos
#     return 1/((A[1]^2+A[2]^2)*(B[1]^2+B[2]^2)*(C[1]^2+C[2]^2))*((A[2]*B[1]-A[1]*B[2])*(C[1]^2+C[2]^2)+(A[1]*C[2]-A[2]*C[1])*(B[1]^2+B[2]^2)+(B[2]*C[1]-B[1]*C[2])*(A[1]^2+A[2]^2)) <=0
# end

function settopology!(e::Edge; new_cwo::Union{Edge,Nothing}=nothing,
    new_ccwo::Union{Edge,Nothing}=nothing, new_cwd::Union{Edge,Nothing}=nothing,
    new_ccwd::Union{Edge,Nothing}=nothing, new_fr::Union{Face,Nothing}=nothing,
    new_fl::Union{Face,Nothing}=nothing)::Nothing
    if !isnothing(new_cwo) e.cwo = new_cwo end
    if !isnothing(new_ccwo) e.ccwo = new_ccwo end
    if !isnothing(new_cwd) e.cwd = new_cwd end
    if !isnothing(new_ccwd) e.ccwd = new_ccwd end
    if !isnothing(new_fr) e.fr = new_fr end
    if !isnothing(new_fl) e.fl = new_fl end
    return
end
settopology!(e::Edge)::Nothing = settopology!(e, new_cwo=e, new_ccwo=e, new_cwd=e, new_ccwd=e)

function getstrut(v::Vertex)::Edge #TODO: add sanity
    @assert !v.original "v is not a dummy vertex"
    edge = v.edge
    while true
        if isstrut(edge)
            return edge
        end
        edge = ccw(edge, v)
    end
    return
end

function findframe(f::Face; dir::Symbol=:ccw)::Edge
	edge = f.edge
	while true
		if isframe(edge)
			return edge
		end
		edge = dir==:ccw ? ccwface(f, edge) : cwface(f, edge)
	end
	return
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
    return all(.!endfield(e, :original))
end

function isstrut(e::Edge)::Bool
    a = endfield(e, :original)
    return xor(a[1], a[2])
end

#Magic
function squeezeedge!(v::Vertex, e::Edge; previous::Union{Edge,Nothing}=nothing, next::Union{Edge,Nothing}=nothing)::Nothing
    @assert v in endpoints(e)
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
    else
        @assert ccw(previous, v) == next
    end
    while true
        if !next_set
            next = ccw(previous, v)
        end
        if previous_set || next_set || sorted_ccw(previous, e, next)
            fr = v==previous.orig ? previous.fl : previous.fr
            fl = v==next.orig ? next.fr : next.fl
            if v == e.orig
                settopology!(e, new_cwo=previous, new_ccwo=next, new_fr=fr, new_fl=fl)
            else
                settopology!(e, new_cwd=previous, new_ccwd=next, new_fr=fl, new_fl=fr)
            end
            cwset!(next, v, e)
            ccwset!(previous, v, e)
            return
        else
            previous = next
        end
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

function endfield(e::Edge, f::Symbol)::Any
    return getfield.(endpoints(e), f)
end

edgelength(e::Edge) = norm(mean(endfield(e, :pos)))::Float64

edgecenter(e::Edge) = mean(endfield(e, :pos))::Array

function joinvertices!(D::DCEL, u::Vertex, v::Vertex, split_face::Bool=true; po::Union{Edge,Nothing}=nothing, pd::Union{Edge,Nothing}=nothing, no::Union{Edge,Nothing}=nothing, nd::Union{Edge,Nothing}=nothing)::Edge
    @assert u != v "Cannot join a vertex to itself"
    new_edge = createedge!(D, u, v)
    squeezeedge!(u, new_edge, previous=po, next=no)
    squeezeedge!(v, new_edge, previous=pd, next=nd)
    disconnected = isnothing(u.edge) || isnothing(v.edge)
    u.edge = v.edge = new_edge
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

function faceverticescw(f::Face)
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

function cwface(f::Face, edge::Union{Edge,Nothing}=nothing)::Edge
    if isnothing(edge)
        return f.edge
    end
    if edge.fr == edge.fl
        return edge.ccwo != edge ? edge.ccwo : edge.ccwd
    elseif f == edge.fr
        return edge.ccwd
    elseif f == edge.fl
        return edge.ccwo
    else #sanity
        throw("$f has a topology problem with edge $edge")
    end
    return
end

function faceedgesccw(f::Face)::Array
    out = []
    edge = f.edge
    while true
        push!(out, edge)
        edge = ccwface(f, edge)
        if edge == f.edge
            break
        end
    end
    return out
end

function faceedgescw(f::Face)::Array
    out = []
    edge = f.edge
    while true
        push!(out, edge)
        edge = cwface(f, edge)
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

function firstlastdummy(D::DCEL)
    v1 = findfirst(x -> !getfield(x, :original), D.vertexlist)
    v2 = findlast(x -> !getfield(x, :original), D.vertexlist)
    return isnothing(v1) ? (nothing, nothing) : (D.vertexlist[v1], D.vertexlist[v2])
end

function findinfiniteface(D::DCEL)::Union{Face,Nothing}
    v1, v2 = firstlastdummy(D)
    if isnothing(v1)
        return nothing
    end
    if v1 == v2
        return v1.edge.fr
    else
        v = v1
        edge = v.edge
        while !isframe(edge) #TODO: add loop protection
            edge = ccw(edge, v)
        end
        fr = edge.fr
        fl = edge.fl
        if isframe(ccwface(fr, edge))
            return fr
        elseif isframe(ccwface(fl, edge))
            return fl
        end
    end
    return
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

function resetface(face::Face, new_face::Face)
    for e in faceedgesccw(face)
        if e.fr == face
            e.fr = new_face
        else
            e.fl = new_face
        end
    end
    return
end

function resetfacelist!(a::DCEL, face::Face, new_face::Face)
    resetface(face, new_face)
    filter!(x->x !=face, a.facelist)
    return
end

function getfaceframe(f::Face, edge::Union{Edge,Nothing}=nothing; dir::Symbol=:ccw)
    isnothing(edge) ? edge=f.edge : @assert f==edge.fr || f==edge.fl
    starting_edge = edge
    while true
        if isframe(edge)
            return edge
        end
        edge = dir==:ccw ? ccwface(f,edge) : cwface(f,edge)
        if edge == starting_edge
            return nothing
        end
    end
    return
end

function addray!(D::DCEL, u::Vertex, angle::Number, f::Union{Face,Nothing}=nothing, split_face::Bool=true)
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
            while true
                outeredge = checkedgefaces(edge)
                if !isnothing(outeredge)
                    f = isframe(ccwface(outeredge.fr, outeredge)) ? outeredge.fl : outeredge.fr
                    break
                end
                edge = ccw(edge, u)
            end
        else
            outeredge = getfaceframe(f)
        end
        e1, e2 = splitedge!(D, outeredge, v)
        previous = ccwface(f, e1) == e2 ? e2 : e1
        new_edge = joinvertices!(D, u, v, split_face, pd=previous)
    end
    return new_edge
end


function vertexinface(v::Vertex, f::Face)
    edge = f.edge
    while true
        if v in endpoints(edge)
            return true
        else
            edge = ccwface(f, edge)
        end
        if edge == f.edge
            return false
        end
    end
    return
end

function checkedgefaces(e::Edge)
    a = getfaceframe(e.fl, e)
    b = getfaceframe(e.fr, e)
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

    new_cwd = e.cwd==e ? new_edge : e.cwd
    new_ccwd = e.ccwd==e ?  new_edge : e.ccwd
    settopology!(new_edge, new_cwo=e, new_ccwo=e, new_cwd=new_cwd, new_ccwd=new_ccwd, new_fr=e.fr, new_fl=e.fl)
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

function leftofline(a::Array, b::Array, c::Array)::Bool
    return ((b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1])) > 0
end

function leftofedge(e::Edge, pos::Array)::Bool
    x, y = e.dest.pos - e.orig.pos
    return cosatan(y,x)*(pos[2]-e.orig.pos[2])-sinatan(y,x)*(pos[1]-e.orig.pos[1]) > 0
end

leftofedge(e::Edge, v::Vertex)::Bool = leftofedge(e, v.pos)

function midpoint(f1::Face, f2::Face)::Array
    return mean([f1.site, f2.site])
end

function perpangle(start::Array, finish::Array)::Float64
    return atan(finish[2]-start[2], finish[1]-start[1]) + pi/2
end
perpangle(f1::Face, f2::Face)::Float64 = perpangle(f1.site, f2.site)

function unstickedge!(e::Edge, v::Vertex)
    if v == e.orig
        ccwset!(e.cwo, e.orig, e.ccwo)
        cwset!(e.ccwo, e.orig, e.cwo)
    elseif v == e.dest
        ccwset!(e.cwd, e.dest, e.ccwd)
        cwset!(e.ccwd, e.dest, e.cwd)
    else
        throw("Given vertex is not in endpoints of edge")
    end
    return
end

function hideedge!(e::Edge)::Nothing
    e.dead = true
    if e.fr.edge == e
        e.fr.edge = e.ccwo
    end
    if e.fl.edge == e
        e.fl.edge = e.ccwo
    end
    unstickedge!(e, e.orig)
    unstickedge!(e, e.dest)
    return
end

function bisectorangles(u::Array, v::Array, w::Array)
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
    addray!(D, m, angle)
    addray!(D, m, angle+pi)
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

function findextrema2(D::DCEL)
    positions = getfield.(D.facelist, :site)
    index_high = 1
    index_low = 1
    for i in 1:length(positions)
        if !isnothing(positions[i])
            if positions[i][2] > positions[index_high][2]
                index_high = i
            end
            if positions[i][2] < positions[index_low][2]
                index_low = i
            end
        end
    end
    return D.facelist[index_high], D.facelist[index_low]
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
    f = max ? findmax : findmin
    index = f(getfield.(getframeface.(edges), :site)[1])[2]
    return edges[index]
end

function findsupport(left::DCEL, right::DCEL, top::Bool)::Tuple#TODO: Optimize
    fl = findinfiniteface(left)
    fr = findinfiniteface(right)
    edge_l = findside(left, true)
    edge_r = findside(right, false)
    f_l = top ? cwface : ccwface
    f_r = top ? ccwface : cwface
    next_l = f_l(fl, edge_l)
    next_r = f_r(fr, edge_r)
    cond = top ? !leftofline : leftofline
    while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_l)) || cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_r))
        while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_l))
            edge_l = next_l
            next_l = f_l(fl, edge_l)
        end
        while cond(getframesite(edge_r), getframesite(edge_l), getframesite(next_r))
            edge_r = next_r
            next_r = f_r(fr, edge_r)
        end
    end
    return getframeface(edge_l), getframeface(edge_r)
end

function findextrema(left::DCEL, right::DCEL)::Tuple
    return (findsupport(left, right, true)..., findsupport(left, right, false)...)
end

function bisector(f1::Face, f2::Face)
    u = f1.site
    v = f2.site
    D = DCEL()
    m = createvertex!(D, mean([u,v]))
    angle = atan(v[2]-u[2], v[1]-u[1]) + pi/2
    addray!(D, m, angle)
    addray!(D, m, angle+pi)
    return D
end

function faceintersection(f::Face)
    edge = f.edge
    while true
        edge = ccwface(f, edge)
        if edge == f.edge
            return false
        end
    end
end

function edgerayintersection(edge::Edge, q::Array, angle::Float64, infinite::Bool=false)::Array
    p = edge.orig.pos
    r = (edge.dest.pos-p)/norm(edge.dest.pos-p)
    s = [cos(angle),sin(angle)]
    num = q - p
    den = cross2d(r, s)
    t = cross2d(num, s)/den
    u = cross2d(num, r)/den
    if isstrut(edge)
        return (u>=0 || infinite) && 0<=t ? p+t*r : [NaN,NaN]
    else
        return (u>=0 || infinite) && 0<=t<=1 ? p+t*r : [NaN,NaN]
    end
    return
end

function rayrayintersection(p::Array, q::Array, v::Array, w::Array)::Array
    x1, y1 = v - p
    x2, y2 = w - q
    r = [cosatan(y1,x1), sinatan(y1,x1)]
    s = [cosatan(y2,x2), sinatan(y2,x2)]
    num = q - p
    den = cross2d(r, s)
    if den==0
        return [NaN,NaN]
    end
    t = cross2d(num, s)/den
    u = cross2d(num, r)/den
    return u>=0 && t>=0 ? p+t*r : [NaN,NaN]
end
function rayrayintersection(p::Vertex, q::Vertex, v::Vertex, w::Vertex)::Array
    return rayrayintersection(p.pos, q.pos, v.pos, w.pos)
end

function rayrayintersection(a::Edge, b::Edge)
    @assert isstrut(a) && isstrut(b)
    return rayrayintersection(a.orig, b.orig, a.dest, b.dest)
end

function rayrayintersection(p::Array, q::Array, a1::Float64, a2::Float64)::Array
    r = [cos(a1),sin(a1)]
    s = [cos(a2),sin(a2)]
    num = q - p
    den = cross2d(r, s)
    if den==0
        return [NaN,NaN]
    end
    t = cross2d(num, s)/den
    u = cross2d(num, r)/den
    return u>=0 && t>=0 ? p+t*r : [NaN,NaN]
end

function facerayintersection(f::Face, start::Array, angle::Number; dir::Symbol=:ccw, infinite::Bool=false, ignore::Array=[], io=stdout)::Tuple
    starter_edge = f.edge
    edge = starter_edge
    write(io, "CHECKING INTERSECTION OF FACE $(f)\n")
    while true
        write(io, "CHECKING EDGE $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        write(io, "WITH LEFT FACE: $(edge.fl) AND RIGHT FACE: $(edge.fr)\n")
        if !isframe(edge) && !(edge in ignore)
            inter = edgerayintersection(edge, start, angle, infinite)
            if !isnan(inter[1])
                write(io, "INTERSECTED $(edge.id) ($(edge.orig.pos), $(edge.dest.pos)), at $inter\n")
                return edge, inter
            end
        end
        next_edge = dir==:cw ? cwface(f,edge) : ccwface(f,edge)
        if next_edge == starter_edge || next_edge == edge
            write(io, "FOUND NO INTERSECTION\n")
            return nothing, [NaN, NaN]
        end
        edge = next_edge
    end
    return
end

function oppositeface(f::Face, e::Edge; io=stdout)::Face
    write(io, "CHECKING OPPOSITE FACE OF: $(f)\n")
    write(io, "WITH MIRROR EDGE: $(e.id) ($(e.orig.pos), $(e.dest.pos))\n")
    write(io, "MIRROR RIGHT FACE: $(e.fr)\n")
    write(io, "MIRROR LEFT FACE: $(e.fl)\n")
    if e.fr.site == f.site
        write(io, "FOUND OPPOSITE: $(e.fl)\n")
        return e.fl
    elseif e.fl.site == f.site
        write(io, "FOUND OPPOSITE: $(e.fr)\n")
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

end
