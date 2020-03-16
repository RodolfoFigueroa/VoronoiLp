##
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

##
mutable struct Handler
    top_left_face::Union{Face,Nothing}
    top_right_face::Union{Face,Nothing}
    current_left_face::Union{Face,Nothing}
    current_right_face::Union{Face,Nothing}
    current_left_vertex::Union{Vertex,Nothing}
    current_right_vertex::Union{Vertex,Nothing}
    current_vertex::Union{Vertex,Nothing}
    current_joint::Union{Edge,Nothing}
    side::Union{Symbol,Nothing}
    ignore::Array
end

##
Handler(lf::Face, rf::Face, lv::Vertex, rv::Vertex) = Handler(lf, rf, lf, rf, lv, rv, nothing, nothing, nothing, [])
##
global vor_sorted = false
function voronoi(start_points::Array)
    if !vor_sorted
        points = sort(start_points, by=x->x[1])
        vor_sorted = true
    else
        points = start_points
    end
    if length(points) == 2
        return voronoitwopoints(points[1], points[2])
    elseif length(points) == 3
        return voronoithreepoints(points[1], points[2], points[3])
    else
        split = Int(floor(length(points)/2))
        points_left = points[1:split]
        points_right = points[split+1:end]
        vor_left = voronoi(points_left)
        vor_right = voronoi(points_right)
        vor = mergevoronoi(vor_left, vor_right)
        return vor
    end
    return
end

function diagramintersection(D::DCEL, il::Array, ir::Array, el::Union{Edge,Nothing}, er::Union{Edge,Nothing}, reference::Union{Vertex,Nothing}=nothing; io=stdout)
    if ir[2] > il[2] || isnan(il[2])
        write(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)++\n")
        starter_right = true
        split = createvertex!(D, ir)
        s = splitedge!(D, er, split)
    elseif il[2] > ir[2] || isnan(ir[2])
        write(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT $(il)++\n")
        starter_right = false
        split = createvertex!(D, il)
        s = splitedge!(D, el, split)
    else #sanity
        throw("")
    end
    return starter_right, split, s
end

function highestintersection(D::DCEL, left::Face, right::Face, start_vertex::Union{Vertex,Nothing}=nothing; io=stdout)
    if isnothing(start_vertex)
        start = midpoint(left, right)
        infinite = true
    else
        start = start_vertex.pos
        infinite = false
    end
    angle = perpangle(right, left)
    write(io, "CHECKING LEFT DIAGRAM...\n")
    el, il = facerayintersection(left, start, angle, dir=:ccw, infinite=infinite, io=io)
    write(io, "CHECKING RIGHT DIAGRAM...\n")
    er, ir = facerayintersection(right, start, angle, dir=:cw, infinite=infinite, io=io)
    dr = distance(start, ir)
    dl = distance(start, il)
    if all(isnan.(il)) || ir[2]>il[2]
        write(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)++\n")
        starter_right = :right
        split = createvertex!(D, ir)
        s = splitedge!(D, er, split)
    elseif all(isnan.(ir)) || il[2]>ir[2]
        write(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT $(il)++\n")
        starter_right = :left
        split = createvertex!(D, il)
        s = splitedge!(D, el, split)
    else #sanity
        throw("")
    end
    return starter_right, split, s
end

function highestintersection(D::DCEL, handler::Handler; io=stdout)
    return highestintersection(D, handler.current_left_face, handler.current_right_face, handler.current_vertex, io=io)
end

function killface!(face::Face, edge::Edge, dir::Symbol, stop_vertex::Union{Vertex,Nothing}=nothing; io=stdout)
    write(io, "KILLING FACE: $face\n")
    while true
        next_edge = dir==:cw ? cwface(face, edge) : ccwface(face, edge)
        edge.dead = true
        v = commonvertex(edge, next_edge)
        write(io, "KILLED EDGE: e$(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        if !v.original
            write(io, "KILLED VERTEX: $v\n")
            v.dead = true
        end
        if stop_vertex in endpoints(edge)
            return
        end
        write(io, "KILLED VERTEX: $v\n")
        v.dead = true
        edge = next_edge
    end
    return
end

function updatehandler!(handler::Handler, edge::Edge, new_vertex::Vertex, side::Symbol; io=stdout)
    if side == :right
        killface!(handler.current_right_face, edge, :cw, handler.current_right_vertex, io=io)
        handler.current_right_face = oppositeface(handler.current_right_face, edge, io=io)
        handler.current_right_vertex = new_vertex
    elseif side == :left
        killface!(handler.current_left_face, edge, :ccw, handler.current_left_vertex, io=io)
        handler.current_left_face = oppositeface(handler.current_left_face, edge, io=io)
        handler.current_left_vertex = new_vertex
    else
        throw("")
    end
    handler.current_vertex = new_vertex
    handler.side = side
    handler.current_joint = new_vertex.edge
    return
end

function foo(D::DCEL, left_face::Face, right_face::Face, left_vertex::Vertex, right_vertex::Vertex; io=stdout)
    handler = Handler(left_face, right_face, left_vertex, right_vertex)
    right_first, split, s = highestintersection(D, handler, io=io)
    if right_first == :right
        edge = cwface(right_face, s[1]) == s[2] ? s[2] : s[1]
    else
        edge = ccwface(left_face, s[1]) == s[2] ? s[2] : s[1]
    end
    updatehandler!(handler, edge, split, right_first, io=io)
    return handler
end

function foobar(D::DCEL, handler::Handler; io=stdout)
    right_first, split, s = highestintersection(D, handler, io=io)
    joint = joinvertices!(D, handler.current_vertex, split, false)
    settopology!(joint, new_fr=handler.current_left_face, new_fl=handler.current_right_face)
    write(io, "CREATED JOINT $(joint.id): ($(joint.orig.pos), $(joint.dest.pos))\n")
    handler.current_left_face.edge = handler.current_right_face.edge = joint
    edge = getfield(joint, right_first==:right ? :ccwd : :cwd)
    if handler.side == :right
        joint.cwo = handler.current_joint
        ccwset!(handler.current_joint, handler.current_vertex, joint)
    else
        joint.ccwo = handler.current_joint
        cwset!(handler.current_joint, handler.current_vertex, joint)
    end
    updatehandler!(handler, edge, split, right_first, io=io)
    return
end

function openface(f::Face, side::Symbol, bearing::Symbol; io=stdout)::Tuple
    e = findframe(f)
    write(io, "OPENED: $e\n")
    e1 = endpoints(e)
    if side == :left
        next_edge = bearing==:top ? ccwface(f, e) : cwface(f, e)
        v = commonvertex(e, next_edge)
        if v == e.orig
            return e, :dest, :orig
        else
            return e, :orig, :dest
        end
    elseif side == :right
        next_edge = bearing==:top ? cwface(f, e) : ccwface(f, e)
        v = commonvertex(e, next_edge)
        if v == e.orig
            return e, :dest, :orig
        else
            return e, :orig, :dest
        end
    else
        throw("")
    end
    return
end

function createray!(D::DCEL, u::Vertex, left::Face, right::Face)
    angle = perpangle(left, right)
    v = createdummyvertex!(D, u, angle)
    joint = joinvertices!(D, u, v, false)
    joint.fl = left
    joint.fr = right
    left.edge = right.edge = joint
    v.edge = joint
    return v
end

function helper(u::Vertex, joint::Edge, t::Tuple, side::Symbol)
    v = getfield(t[1],t[3])
    if side == :left
        temp = t[3] == :orig ? t[1].cwo : t[1].cwd
        j = joinvertices!(D, u, v, false, pd=temp, no=joint)
    else
        temp = t[3] == :orig ? t[1].ccwo : t[1].ccwd
        j = joinvertices!(D, u, v, false, nd=temp, po=joint)
    end
    hideedge!(t[1])
    return j
end

function mergevoronoi(left::DCEL, right::DCEL, logfile::IOStream)
    io = isnothing(logfile) ? stdout : logfile
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    push!(D.facelist, fi)
    hl, hr, ll, lr = findextrema(left, right)
    write(io, "TOP FACES: ($hl, $hr)\n")
    write(io, "BOT FACES: ($ll, $lr)\n")

    tl, tr, bl, br = openface.([hl, hr, ll, lr], [:left, :right, :left, :right], [:top, :top, :bot, :bot], io=io)
    starter_left = getfield(tl[1], tl[2])
    starter_right = getfield(tr[1], tr[2])
    global handler = foo(D, hl, hr, starter_left, starter_right, io=io)

    top_vertex = handler.current_vertex
    top = createray!(D, top_vertex, hl, hr)
    joint = top.edge
    handler.current_joint = joint
    j = helper(top, joint, tl, :left)
    j.fr = fi
    fi.edge = j
    j = helper(top, joint, tr, :right)
    j.fl = fi

    bottom_vertex = nothing
    while true
        if handler.current_left_face == ll && handler.current_right_face == lr
            bottom_vertex = handler.current_vertex
            break
        end
        write(io, "\n%%%TOP METHOD%%%\n")
        foobar(D, handler, io=io)
    end
    bot = createray!(D, bottom_vertex, lr, ll)
    if handler.side == :right
        bot.edge.cwo = handler.current_joint
        handler.current_joint.ccwd = bot.edge
    else
        bot.edge.ccwo = handler.current_joint
        handler.current_joint.cwd = bot.edge

    end
    joint = bot.edge
    j = helper(bot, joint, bl, :right)
    j.fl = fi
    j = helper(bot, joint, br, :left)
    j.fr = fi
    # filter!(x->!getfield(x, :dead), D.edgelist)
    # filter!(x->!getfield(x, :dead), D.vertexlist)
    # fixids!(D)
    return D
end
##
# p = [[2,2],[1,1],[3,5]]
# p = [[0.59,0.58],[0.32,0.92],[0.45,0.71]]
p = [[0.87, 0.60], [0.08,0.61], [0.73,0.63]]
q = [[3,3],[5,4],[6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
fixids!(a)
fixids!(b)
checkdcel(a)
checkdcel(b)
# p1 = plotdcel(a, faces=true)
# p2 = plotdcel(b, faces=true)
# plot(p1, p2)

file = open("log.txt", "w")
test = mergevoronoi(a,b,file)
close(file)
# plotdcel(test)

##
plotdcel(D, ratio=:equal, sites=false, dead_edges=false, dead_vertices=false)
# checkdcel(test)

##
# splittest = createvertex!(a, [2.65, 2.48])
# splitedge!(a, a.edgelist[2], splittest)
# addray!(a, splittest, deg2rad(-56)+pi, a.facelist[4], false)
# plotdcel(a)
