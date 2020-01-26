##
using Revise
include("./vector.jl")
includet("./struct.jl")
using .DVector
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

function killface!(edge::Edge, joint::Edge, right::Bool; io=stdout)
    while true
        next_edge = right ? cwface(joint.fl, edge) : ccwface(joint.fr, edge)
        hideedge(edge)
        v = commonvertex(edge, next_edge)
        v.dead = true
        write(io, "KILLED VERTEX: $(v.id)")
        write(io, "KILLED (1): $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))")
        if edge == next_edge #edge.dead #|| current_left_vertex in endpoints(edge)
            break
        end
        edge = next_edge
    end
    return
end

function firstray!(split::Vertex, angle::Number, hr::Face, hl::Face, right::Bool; io=stdout)
    upper_ray = addray!(D, split, angle+pi, right ? hr : hl, false)
    write(io, "RAY HAS AN ORIGIN $(upper_ray.orig.pos)\n")
    starter_edge = right ? ccw(upper_ray, split) : cw(upper_ray, split)
    edge = starter_edge
    while true
        hideedge(edge)
        write(io, "KILLED (2): $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        if isboundaryedge(edge)
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


function diagramintersection(ir::Array, il::Array, er::Union{Edge,Nothing}, el::Union{Edge,Nothing}; io=stdout)
    if ir[2] > il[2]
        write(io, "\n++INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)++\n")
        starter_right = true
        split = createvertex!(D, ir)
        s1, s2 = splitedge!(D, er, split)
    elseif il[2] > ir[2]
        write(io, "\n++INTERSECTED LEFT DIAGRAM FIRST AT $(il)++\n")
        starter_right = false
        split = createvertex!(D, il)
        s1, s2 = splitedge!(D, el, split)
    else #sanity
        throw("")
    end
    return starter_right, split, s1, s2
end

function weldedges(joint::Edge, lr::Face, ll::Face, right::Bool; io=stdout)
    write(io, "WELDING...\n")
    weld_edge = right ? joint.cwo : joint.ccwo
    next_weld_edge = Edge
    while true
        next_weld_edge = right ? ccwface(ll, weld_edge) : cwface(lr, weld_edge)
        if isboundaryedge(next_weld_edge)
            common_weld_vertex = commonvertex(weld_edge, next_weld_edge)
            weld_vertex = common_weld_vertex==next_weld_edge.orig ? next_weld_edge.dest : next_weld_edge.orig
            return next_weld_edge, weld_vertex
        end
        weld_edge = next_weld_edge
    end
    return
end

function upperweld(edge::Edge, upper_ray::Edge, hr::Face, hl::Face, right::Bool; io=stdout)
    while true
        hideedge(edge)
        write(io, "KILLED (3): $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))\n")
        next_edge = right ? cwface(hr, edge) : ccwface(hl, edge)
        commonvertex(edge, next_edge).dead = true
        if isboundaryedge(next_edge)
            common_vertex = commonvertex(edge, next_edge)
            unstickedge!(next_edge, common_vertex)
            if common_vertex == next_edge.dest
                next_edge.dest = upper_ray.dest
                right ? next_edge.fl = hr : next_edge.fr = hl
            elseif common_vertex == next_edge.orig
                next_edge.orig = upper_ray.dest
                right ? next_edge.fr = hr : next_edge.fl = hl
            end
            if right
                squeezeedge!(upper_ray.dest, next_edge, previous=upper_ray)
            else
                squeezeedge!(upper_ray.dest, next_edge, next=upper_ray)
            end
            break
        end
        edge = next_edge
    end
    return
end

#When I wrote this, only God and I understood what I was doing.
#Now, only God knows.
function mergevoronoi(left::DCEL, right::DCEL, logging::Bool=false)
    global io = logging ? open("log.txt", "w") : stdout
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    push!(D.facelist, fi)
    hl, ll = findextrema(left)
    hr, lr = findextrema(right)
    angle = perpangle(hr, hl)
    write(io, "==INITIAL RAY HAS AN ANGLE: $(rad2deg(angle))==\n")
    write(io, "CHECKING LEFT DIAGRAM\n")
    el, il = facerayintersection(hl, midpoint(hr, hl), angle, false, io=io)
    write(io, "CHECKING LEFT DIAGRAM\n")
    er, ir = facerayintersection(hr, midpoint(hr, hl), angle, true, io=io)
    starter_right, split, s1, s2 = diagramintersection(ir, il, er, el, io=io)

    current_right_face = current_left_face = Face
    current_right_vertex = current_left_vertex = nothing
    upper_ray = Edge
    if starter_right
        upper_ray = firstray!(split, angle, hr, hl, true, io=io)
        current_right_face = oppositeface(hr, upper_ray.cwo, io=io)
        current_left_face = hl
        current_right_vertex = split
    else
        upper_ray = firstray!(split, angle, hr, hl, false, io=io)
        current_right_face = hr
        current_left_face = oppositeface(hl, upper_ray.ccwo, io=io)
        current_left_vertex = split
    end
    current_vertex = split
    weld_set = false
    weld_vertex = Vertex
    next_weld_edge = Edge
    finisher_right = Bool
    while true
        angle = perpangle(current_right_face, current_left_face)
        write(io, "\n==NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $(current_vertex.pos)==\n")
        if current_left_face == ll && current_right_face == lr
            lower_ray = addray!(D, current_vertex, angle, finisher_right ? lr : ll, false)
            settopology!(lower_ray, new_fr=ll, new_fl=lr)
            hideedge(finisher_right ? lower_ray.ccwd : lower_ray.cwd)
            write(io, "KILLED (4): $(lower_ray.ccwd)\n")
            if next_weld_edge.orig == weld_vertex
                next_weld_edge.orig = lower_ray.dest
            elseif next_weld_edge.dest == weld_vertex
                next_weld_edge.dest = lower_ray.dest
            else #sanity
                throw("")
            end
            if finisher_right
                squeezeedge!(lower_ray.dest, next_weld_edge, previous=lower_ray)
            else
                squeezeedge!(lower_ray.dest, next_weld_edge, next=lower_ray)
            end
            filter!(x -> !getfield(x, :dead), D.edgelist)
            filter!(x -> !getfield(x, :dead), D.vertexlist)
            fixids!(D)
            if logging
                close(io)
            end
            return D
        end
        if current_right_face != lr
            er, ir = facerayintersection(current_right_face, current_vertex.pos, angle, true, [s1,s2], io=io)
        else
            ir = [-Inf, -Inf]
        end
        if current_left_face != ll
            el, il = facerayintersection(current_left_face, current_vertex.pos, angle, false, [s1,s2], io=io)
        else
            il = [-Inf, -Inf]
        end
        current_edge = Edge
        current_split = Vertex
        right_first = Bool
        current_edge = ir[2]>il[2] ? er : el
        right_first, current_split, s1, s2 = diagramintersection(ir, il, er, el, io=io)
        joint = joinvertices!(D, current_vertex, current_split, false)
        settopology!(joint, new_fr=current_left_face, new_fl=current_right_face)
        write(io, "CREATED JOINT $(joint.id): ($(joint.orig.pos), $(joint.dest.pos)) AND ASSIGNED FL:$current_right_face, FR:$current_left_face\n")
        current_left_face.edge = current_right_face.edge = joint
        edge = right_first ? joint.ccwd : joint.cwd
        write(io, "CURRENTLY ON EDGE ($(edge.orig.pos), $(edge.dest.pos))\n")
        if right_first
            if !weld_set && current_left_face == ll
                next_weld_edge, weld_vertex = weldedges(joint, lr, ll, true, io=io)
                finisher_right = true
                weld_set = true
                killface!(edge, joint, true, io=io)
            elseif current_right_face == hr
                upperweld(edge, upper_ray, hr, hl, true, io=io)
            else
                killface!(edge, joint, true, io=io)
            end
            current_right_face = oppositeface(current_right_face, current_edge, io=io)
        else
            if !weld_set && current_right_face == lr
                next_weld_edge, weld_vertex = weldedges(joint, lr, ll, false, io=io)
                finisher_right = false
                weld_set = true
                killface!(edge, joint, false, io=io)
            elseif current_left_face == hl
                upperweld(edge, upper_ray, hr, hl, false, io=io)
            else
                killface!(edge, joint, false, io=io)
            end
            current_left_face = oppositeface(current_left_face, current_edge, io=io)
        end
        current_vertex = current_split
    end
    return
end


##
# p = [[2,2],[1,1],[3,5]]
# p = [rand(2), rand(2), rand(2)]
q = [[3,3],[5,4],[6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
fixids!(a)
fixids!(b)
checkdcel(a)
checkdcel(b)
test = mergevoronoi(a, b, true)
checkdcel(test)

##
plotdcel(test, ratio=:equal, sites=true, dead_edges=true, dead_vertices=false)
# checkdcel(test)
