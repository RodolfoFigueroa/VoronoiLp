module DStruct
include("./vector.jl")
import Base.show
using Plots, Statistics, .DVector
export Vertex, Face, Edge, DCEL, show, createvertex!, createedge!, createface!, createdummyvertex!, fixids!, checkdcel, plotdcel, cw, ccw, cwset!, ccwset!, add_join_vertex!, joinvertices!, addray!, splitedge!, deleteedge!, endpoints, midpoint, isboundaryedge, isstrutedge, cwface, ccwface, commonvertex, resetfacelist!, leftofedge, perpangle, settopology!, squeezeedge!

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
    cwo::Union{Edge,Nothing}
    ccwo::Union{Edge,Nothing}
    cwd::Union{Edge,Nothing}
    ccwd::Union{Edge,Nothing}
    fr::Union{Face,Nothing}
    fl::Union{Face,Nothing}
    dead::Bool
end
Edge(id, orig, dest) = Edge(id, orig, dest, nothing, nothing, nothing, nothing, nothing, nothing, false)

mutable struct DCEL
    edgelist::Array
    vertexlist::Array
    facelist::Array
    dummylist::Array
end
DCEL() = DCEL([], [], [], [])

function printnothing(e::Any)::String
    p = String
    if typeof(e) == Edge
        p = "e"
    elseif typeof(e) == Face
        p = "f"
    elseif typeof(e) == Vertex
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
    push!(D.dummylist, v)
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

function plotdcel(D::DCEL; faces::Bool=false, dead::Bool=false, bounds::Bool=true, ratio::Any=:equal, sites::Bool=false)
    p = plot(leg=false)
    for v in D.vertexlist
        color = v.original ? :black : :red
        scatter!(p, [v.pos[1]], [v.pos[2]], series_annotations=[Plots.text("\nv$(v.id)"), :bottom], color=color, aspect_ratio=ratio)
    end
    for e in D.edgelist
        bound = isboundaryedge(e)
        if e.dead
            color = :gray
        else
            color = bound ? :red : :black
        end
        ave = edgecenter(e)
        if (!e.dead || dead) && (!bound || bounds)
            plot!(p, [e.orig.pos[1],e.dest.pos[1]], [e.orig.pos[2],e.dest.pos[2]], line=:arrow, annotations = (ave[1], ave[2], "e$(e.id)"), linecolor=color, aspect_ratio=ratio)
        end
    end
    for f in D.facelist
        if faces
            vertices = faceverticesccw(f)
            x,y = mean(getfield.(vertices, :pos))
            annotate!([x],[y], "f$(f.id)")
        end
        if sites
            if !isnothing(f.site)
                scatter!([f.site[1]], [f.site[2]])
            end
        end
    end
    return p
end

function checkdcel(D::DCEL)
    for e in D.edgelist
        @assert cw(ccw(e,e.orig),e.orig) == e "1 $(e)"
        @assert ccw(cw(e,e.orig),e.orig) == e "2 $(e)"
        @assert cw(ccw(e,e.dest),e.dest) == e "3 $(e)"
        @assert ccw(cw(e,e.dest),e.dest) == e "4 $(e)"
        @assert e.orig != e.dest
        if isstrutedge(e)
            @assert e.orig.original && !e.dest.original
        end
    end

    pos_list = getfield.(D.vertexlist, :pos)

    for f in D.facelist
        edges = faceedgesccw(f)
        for e in edges
            @assert f==e.fr || f==e.fl "$(f):$(e)"
        end
        edges = faceedgescw(f)
        for e in edges
            @assert f==e.fr || f==e.fl "$(f):$(e)"
        end

    end
    return
end

cw(e::Edge, v::Vertex)::Edge = v==e.orig ? e.cwo : e.cwd

ccw(e::Edge, v::Vertex)::Edge = v==e.orig ? e.ccwo : e.ccwd

cwset!(e::Edge, v::Vertex, new_edge::Edge) = v==e.orig ? e.cwo=new_edge : e.cwd=new_edge

ccwset!(e::Edge, v::Vertex, new_edge::Edge) = v==e.orig ? e.ccwo=new_edge : e.ccwd=new_edge

endpoints(u::Edge)::Set = Set([u.orig, u.dest])

edgeangle(u::Edge)::Float64 = atan(u.dest.pos[2]-u.orig.pos[2], u.dest.pos[1]-u.orig.pos[1])

pivotangle(u::Edge, p::Vertex)::Float64 = edgeangle(u) + (p==u.orig ? 0 : pi)

function commonvertex(u::Edge, edges...)
    f = endpoints(u)
    for e in endpoints.(edges)
        intersect!(f, e)
    end
    return isempty(f) ? nothing : pop!(f)
end

function uncommonvertices(u::Edge, v::Edge, w::Edge)
    s1 = endpoints(u)
    s2 = endpoints(v)
    s3 = endpoints(w)
    s = intersect(s1, s2, s3)
    return pop!(setdiff(s1,s)), pop!(setdiff(s2,s)), pop!(setdiff(s3,s)), pop!(s)
end

function sorted_ccw(u::Edge, v::Edge, w::Edge)::Bool
    pivot = commonvertex(u, v, w)
    a = pivotangle(u, pivot)
    b = pivotangle(v, pivot)
    c = pivotangle(w, pivot)
    return sin(a-b) - sin(a-c) + sin(b-c) <= 0
end

function sorted_ccw2(u::Edge, v::Edge, w::Edge)::Bool
    a, b, c, d = uncommonvertices(u, v, w)
    return (b.pos[1]-a.pos[1])*(b.pos[2]+a.pos[2])+(c.pos[1]-b.pos[1])*(c.pos[2]+b.pos[2])+(d.pos[1]-c.pos[1])*(d.pos[2]+c.pos[2])+(a.pos[1]-d.pos[1])*(a.pos[2]+d.pos[2]) <= 0
end

function sorted_ccw3(u::Edge, v::Edge, w::Edge)::Bool
    p,q,r,z = uncommonvertices(u, v, w)
    A = p.pos-z.pos
    B = q.pos-z.pos
    C = r.pos-z.pos
    return 1/((A[1]^2+A[2]^2)*(B[1]^2+B[2]^2)*(C[1]^2+C[2]^2))*((A[2]*B[1]-A[1]*B[2])*(C[1]^2+C[2]^2)+(A[1]*C[2]-A[2]*C[1])*(B[1]^2+B[2]^2)+(B[2]*C[1]-B[1]*C[2])*(A[1]^2+A[2]^2)) <=0
end

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

function getstrut(v::Vertex)::Edge
    @assert !v.original "v is not a dummy vertex"
    edge = v.edge
    while true
        if any(endfield(edge, :original))
            return edge
        end
        edge = ccw(edge, v)
    end
    return
end

function normalizedummy(v::Vertex)::Array
    @assert !v.original "v is not a dummy vertex"
    e = getstrut(v)
    u, v = e.orig.original ? (e.orig, e.dest) : (e.dest, e.orig)
    return v.pos - u.pos
end

vertexccw(vertices::Array) = pointccw(getfield.(vertices, :pos))

function isboundaryedge(e::Edge)::Bool
    return all(.!endfield(e, :original))
end

function isstrutedge(e::Edge)::Bool
    a = endfield(e, :original)
    return xor(a[1], a[2])
end

#Magic
function squeezeedge!(v::Vertex, e::Edge, previous::Union{Edge,Nothing}=nothing)::Nothing
    @assert v in endpoints(e)
    if isnothing(v.edge)
        return
    end
    previous_set = true
    if isnothing(previous)
        previous_set = false
        previous = v.edge
    end
    while true
        next = ccw(previous, v)
        if previous_set || sorted_ccw(previous, e, next)
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

endfield(e::Edge, f::Symbol) = getfield.(endpoints(e), f)::Any

edgelength(e::Edge) = norm(mean(endfield(e, :pos)))::Float64

edgecenter(e::Edge) = mean(endfield(e, :pos))::Array

function joinvertices!(D::DCEL, u::Vertex, v::Vertex; p1::Union{Edge,Nothing}=nothing, p2::Union{Edge,Nothing}=nothing)::Edge
    @assert u != v "Cannot join a vertex to itself"
    new_edge = createedge!(D, u, v)
    squeezeedge!(u, new_edge, p1)
    squeezeedge!(v, new_edge, p2)
    disconnected = isnothing(u.edge) || isnothing(v.edge)
    u.edge = v.edge = new_edge
    if disconnected
        return new_edge
    end

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
        throw("Looking for $(f) in $(edge.fl) and $(edge.fr)")
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
    else
        throw("$(f) has a topology problem")
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

function findinfiniteface(D::DCEL)::Union{Face,Nothing}
    if isempty(D.dummylist)
        return nothing
    end
    if length(D.dummylist) == 1
        return D.dummylist[1].edge.fr
    else
        v = D.dummylist[1]
        edge = v.edge
        while !isboundaryedge(edge)
            edge = ccw(edge, v)
        end
        global fr = edge.fr
        global fl = edge.fl
        if isboundaryedge(ccwface(fr, edge))
            return fr
        elseif isboundaryedge(ccwface(fl, edge))
            return fl
        end
    end
    return
end

function deleteedge!(D::DCEL, e::Edge)::Nothing
    f = Face("$(length(D.facelist))")
    f.edge = e.ccwo
    resetfacelist!(D, f, e.fl)
    resetfacelist!(D, f, e.fr)
    ccwset!(e.cwo, e.orig, e.ccwo)
    cwset!(e.ccwo, e.orig, e.cwo)
    ccwset!(e.cwd, e.dest, e.ccwd)
    cwset!(e.ccwd, e.dest, e.cwd)
    filter!(x->x != e, D.edgelist)
    push!(D.facelist, f)
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

function resetfacelist!(a::DCEL, new_face::Face, face::Union{Face,Nothing}=nothing)
    if isnothing(face)
        face = findinfiniteface(a)
    end
    resetface(face, new_face)
    filter!(x->x !=face, a.facelist)
    return
end

function addray!(D::DCEL, u::Vertex, angle::Number, f::Union{Face,Nothing}=nothing)
    l = length(D.dummylist)
    v = createdummyvertex!(D, u, angle)
    new_edge = nothing
    if l <=1
        new_edge = joinvertices!(D, u, v)
        if l == 1
            outer = joinvertices!(D, D.dummylist[1], v)
            joinvertices!(D, v, D.dummylist[1], p1=outer, p2=getstrut(D.dummylist[1]))
        end
    else
        outeredge = nothing
        if isnothing(f)
            edge = u.edge
            while true
                outeredge = checkedgefaces(edge)
                if !isnothing(outeredge)
                    f = isboundaryedge(ccwface(outeredge.fr, outeredge)) ? outeredge.fl : outeredge.fr
                    break
                else
                    edge = ccw(edge, u)
                end
            end
        else
            outeredge = checkfacebounds(f)
        end
        e1, e2 = splitedge!(D, outeredge, v)
        previous = ccwface(f, e1) == e2 ? e2 : e1
        new_edge = joinvertices!(D, u, v, p2=previous)
    end
    return new_edge
end

function checkfacebounds(f::Face, edge::Union{Edge,Nothing}=nothing)
    if isnothing(edge)
        edge = f.edge
    else
        @assert f==edge.fr || f==edge.fl
    end
    starting_edge = edge
    while true
        if isboundaryedge(edge)
            return edge
        end
        edge = ccwface(f, edge)
        if edge == starting_edge
            return nothing
        end
    end
    return
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
    a = checkfacebounds(e.fl, e)
    b = checkfacebounds(e.fr, e)
    if !isnothing(a)
        return a
    elseif !isnothing(b)
        return b
    else
        return nothing
    end
end

function splitedge!(D::DCEL, e::Edge, split::Vertex)
    new_edge = Edge("$(e.id)b", split, e.dest)
    push!(D.edgelist, new_edge)

    new_cwd = e.cwd==e ? new_edge : e.cwd
    new_ccwd = e.ccwd==e ?  new_edge : e.ccwd
    settopology!(new_edge, new_cwo=e, new_ccwo=e, new_cwd=new_cwd, new_ccwd=new_ccwd, new_fr=e.fr, new_fl=e.fl)
    split.edge = new_edge
    cwset!(e.ccwd, e.dest, new_edge)
    ccwset!(e.cwd, e.dest, new_edge)
    if isstrutedge(e)
        x, y = e.dest.pos - e.orig.pos
        new_edge.dest.pos = split.pos + [cosatan(y,x), sinatan(y,x)]
    end

    e.id = "$(e.id)a"
    e.dest.edge = new_edge
    e.dest = split
    e.cwd = e.ccwd = new_edge
    return e, new_edge
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

end
