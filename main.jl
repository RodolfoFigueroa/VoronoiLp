##
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

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

function firstray!(D::DCEL, split::Vertex, angle::Number, hr::Face, hl::Face, right::Bool; io=stdout)
    upper_ray = addray!(D, split, angle+pi, right ? hr : hl, false)
    write(io, "RAY HAS AN ORIGIN $(upper_ray.orig.pos)\n")
    edge = right ? ccw(upper_ray, split) : cw(upper_ray, split)
    while true
        hideedge!(edge)
        write(io, "KILLED (2): $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        if isframe(edge)
            break
        end
        next = right ? cwface(upper_ray.fl, edge) : ccwface(upper_ray.fr, edge)
        commonvertex(edge, next).dead = true
        edge = next
    end
    upper_ray.fl = hl
    upper_ray.fr = hr
    return upper_ray
end

function diagramintersection(D::DCEL, ir::Array, il::Array, er::Union{Edge,Nothing}, el::Union{Edge,Nothing}, side::Symbol; io=stdout)
    if (side==:top && (ir[2] > il[2] || isnan(il[2]))) || (side==:bot && (ir[2] < il[2] || isnan(il[2])))
        write(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)++\n")
        starter_right = true
        split = createvertex!(D, ir)
        s = splitedge!(D, er, split)
    elseif (side==:top && (il[2] > ir[2] || isnan(ir[2]))) || (side==:bot && (il[2] < ir[2] || isnan(ir[2])))
        write(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT $(il)++\n")
        starter_right = false
        split = createvertex!(D, il)
        s = splitedge!(D, el, split)
    else #sanity
        throw("")
    end
    return starter_right, split, s
end

function openface(f::Face)
    hideedge!(findframe(f))
    return
end

function killface!(face::Face, edge::Edge, dir::Symbol, stop_at_frame::Bool=false; io=stdout)
    write(io, "KILLING FACE: $face\n")
    while true
        next_edge = dir==:cw ? cwface(face, edge) : ccwface(face, edge)
        hideedge!(edge)
        v = commonvertex(edge, next_edge)
        v.dead = true
        write(io, "KILLED VERTEX: $v\n")
        write(io, "KILLED EDGE: e$(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        if stop_at_frame
            if isframe(edge)
                v.dead = false
                return v
            end
        else
            if edge == next_edge
                return
            end
        end
        edge = next_edge
    end
    return
end

function foo(D::DCEL, left::Face, right::Face, side::Symbol; io=stdout)
    angle = side==:top ? perpangle(right, left) : perpangle(left, right)
    m = midpoint(right, left)
    write(io, "==INITIAL RAY HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $m==\n")
    write(io, "CHECKING LEFT DIAGRAM...\n")
    el, il = facerayintersection(left, m, angle, dir=:ccw, infinite=true, io=io)
    write(io, "CHECKING RIGHT DIAGRAM...\n")
    er, ir = facerayintersection(right, m, angle, dir=:cw, infinite=true, io=io)
    right_first, split, s = diagramintersection(D, ir, il, er, el, side, io=io)
    top_vertex = current_vertex = split
    current_right_face = current_left_face = Face
    top_left_vertex = top_right_vertex = bottom_right_vertex = bottom_left_vertex = Vertex
    if right_first
        current_right_face = oppositeface(right, s[1], io=io)
        current_left_face = left
        edge = cwface(right, s[1]) == s[2] ? s[2] : s[1]
        top_right_vertex = killface!(right, edge, :cw, true, io=io)
    else
        current_right_face = right
        current_left_face = oppositeface(left, s[1], io=io)
        edge = ccwface(left, s[1]) == s[2] ? s[2] : s[1]
        top_left_vertex = killface!(left, edge, :ccw, true, io=io)
    end
    return current_left_face, current_right_face, current_vertex
end

function foobar(D::DCEL, left::Face, right::Face, top_left::Face, top_right::Face, current_vertex::Vertex, s::Array, side::Symbol; io=stdout)
    angle = side==:top ? perpangle(right, left) : perpangle(left, right)
    write(io, "==NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $(current_vertex.pos)==\n")
    er, ir = facerayintersection(right, current_vertex.pos, angle, dir=:cw, ignore=s, io=io)
    el, il = facerayintersection(left, current_vertex.pos, angle, dir=:ccw, ignore=s, io=io)
    right_first, current_split, s = diagramintersection(D, ir, il, er, el, side, io=io)
    joint = joinvertices!(D, current_vertex, current_split, false)
    if side == :top
        settopology!(joint, new_fr=left, new_fl=right)
    elseif side == :bot
        settopology!(joint, new_fr=right, new_fl=left)
    end
    write(io, "CREATED JOINT $(joint.id): ($(joint.orig.pos), $(joint.dest.pos))\n")
    left.edge = right.edge = joint
    if right_first
        b = right==top_right
        current_edge = joint.ccwd
        temp = killface!(right, current_edge, :cw, b, io=io)
        if b
            top_right_vertex = temp
        end
        right = oppositeface(right, current_edge, io=io)
    else
        b = left==top_left
        current_edge = joint.cwd
        temp = killface!(left, joint.cwd, :ccw, b, io=io)
        if b
            top_left_vertex = temp
        end
        left = oppositeface(left, current_edge, io=io)
    end
    return current_split, left, right, s
end

function mergevoronoi(left::DCEL, right::DCEL, logfile::IOStream)
    io = isnothing(logfile) ? stdout : logfile
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    push!(D.facelist, fi)
    hl, hr, ll, lr = findextrema(left, right)
    write(io, "TOP FACES: ($hl, $hr)\n")
    write(io, "BOTTOM FACES: ($ll, $lr)\n")

    current_top_left_face, current_top_right_face, current_top_vertex = foo(D, hl, hr, :top, io=io)
    current_bottom_left_face, current_bottom_right_face, current_bottom_vertex = foo(D, ll, lr, :bot, io=io)

    s_top = s_bottom = [nothing, nothing]
    bottom_vertex = current_split = Vertex
    while true
        write(io, "\n%%%TOP METHOD%%%\n")
        current_top_vertex, current_top_left_face, current_top_right_face = foobar(D, current_top_left_face, current_top_right_face, hl, hr, current_top_vertex, s_top, :top; io=io)
        if current_top_vertex == current_bottom_vertex
            break
        end
        write(io, "\n%%%BOTTOM METHOD%%%\n")
        current_bottom_vertex, current_bottom_left_face, current_bottom_right_face = foobar(D, current_bottom_left_face, current_bottom_right_face, ll, lr, current_bottom_vertex, s_bottom, :bot; io=io)
        sleep(0.1)
    end
    # println(top_left_vertex)
    # println(top_right_vertex)
    # println(bottom_left_vertex)
    # println(bottom_right_vertex)
    # println(bottom_vertex)
    return D
end
##
p = [[2,2],[1,1],[3,5]]
# p = [[0.59,0.58],[0.32,0.92],[0.45,0.71]]
q = [[3,3],[5,4],[6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
fixids!(a)
fixids!(b)
checkdcel(a)
checkdcel(b)
p1 = plotdcel(a, faces=true)
p2 = plotdcel(b, faces=true)
plot(p1, p2)

file = open("log.txt", "w")
global top, bot = mergevoronoi(a,b,file)
close(file)
# plotdcel(test)

##
plotdcel(test, ratio=:equal, sites=true, dead_edges=true, dead_vertices=false)
# checkdcel(test)

##
# splittest = createvertex!(a, [2.65, 2.48])
# splitedge!(a, a.edgelist[2], splittest)
# addray!(a, splittest, deg2rad(-56)+pi, a.facelist[4], false)
# plotdcel(a)
