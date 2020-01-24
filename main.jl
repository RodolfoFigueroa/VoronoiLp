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

#When I wrote this, only God and I understood what I was doing.
#Now, only God knows.
function mergevoronoi(left::DCEL, right::DCEL)
    counter = 0
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    push!(D.facelist, fi)
    hl, ll = findextrema(left)
    hr, lr = findextrema(right)
    angle = perpangle(hr, hl)
    println("INITIAL RAY HAS AN ANGLE: $(rad2deg(angle))")
    el, il = facerayintersection(hl, midpoint(hr, hl), angle, false)
    er, ir = facerayintersection(hr, midpoint(hr, hl), angle, true)
    split = Vertex
    starter_right = Bool
    s1 = s2 = Edge
    if ir[2] > il[2]
        println("INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)")
        starter_right = true
        split = createvertex!(D, ir)
        s1, s2 = splitedge!(D, er, split)
    else
        println("INTERSECTED LEFT DIAGRAM FIRST AT $(il)")
        starter_right = false
        split = createvertex!(D, il)
        s1, s2 = splitedge!(D, el, split)
    end
    current_right_face = current_left_face = Face
    current_right_vertex = current_left_vertex = nothing
    upper_ray = Edge
    if starter_right
        upper_ray = addray!(D, split, angle+pi, hr, false)
        println("RAY HAS AN ORIGIN $(upper_ray.orig.pos)")
        starter_edge = ccw(upper_ray, split)
        edge = starter_edge
        while true
            # edge.dead = true
            hideedge(edge)
            println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
            if isboundaryedge(edge)
                break
            end
            next = cwface(upper_ray.fl, edge)
            commonvertex(edge, next).dead = true
            edge = next
        end
        current_right_face = oppositeface(upper_ray.fl, starter_edge)
        current_left_face = hl
        current_vertex = split
        upper_ray.fl = hl
        upper_ray.fr = hr
        current_right_vertex = current_vertex
        println("\nCURRENT RIGHT VERTEX IS NOW: $(current_right_vertex.pos)")
    end
    while true
        if current_left_face == ll && current_right_face == lr
            # filter!(x -> !getfield(x, :dead), D.edgelist)
            # filter!(x -> !getfield(x, :dead), D.vertexlist)
            fixids!(D)
            return D
        end
        angle = perpangle(current_right_face, current_left_face)
        println("\n\n==NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $(current_vertex.pos)==")
        if current_right_face != lr
            er, ir = facerayintersection(current_right_face, current_vertex.pos, angle, true, [s1,s2])
            if !isnothing(er)
                println("\n=BISECTOR INTERSECTED RIGHT EDGE: ($(er.orig.pos), $(er.dest.pos) AT: $(ir)=")
            end
        else
            ir = [-Inf, -Inf]
        end
        if current_left_face != ll
            el, il = facerayintersection(current_left_face, current_vertex.pos, angle, false, [s1,s2])
            if !isnothing(el)
                println("\n=BISECTOR INTERSECTED LEFT EDGE ($(el.orig.pos), $(el.dest.pos) AT: $(il)=")
            end
        else
            il = [-Inf, -Inf]
        end
        current_edge = Edge
        current_split = Vertex
        right_first = Bool
        if ir[2] > il[2]
            println("\n+CHOSE RIGHT INTERSECTION+\n")
            right_first = true
            current_split = createvertex!(D, ir)
            current_edge = er
        elseif il[2] > ir[2]
            println("\n+CHOSE LEFT INTERSECTION+\n")
            right_first = false
            current_split = createvertex!(D, il)
            current_edge = el
        else
            throw("")
        end
        println("SPLITTING: ($(current_edge.orig.pos), $(current_edge.dest.pos))")
        s1, s2 = splitedge!(D, current_edge, current_split)
        println("SPLIT RESULTS: $(s1), $(s2)")
        joint = joinvertices!(D, current_vertex, current_split, false)
        joint.fr = current_left_face
        joint.fl = current_right_face
        current_left_face.edge = current_right_face.edge = joint
        if right_first
            edge = ccw(joint, current_split)
            while true
                next_edge = cwface(joint.fl, edge)
                if edge == next_edge || current_right_vertex in endpoints(edge)#edge.dead #|| current_right_vertex in endpoints(edge)
                    current_right_vertex = current_split
                    println("CURRENT RIGHT VERTEX IS: $(current_right_vertex.pos)")
                    break
                end
                # edge.dead = true
                hideedge(edge)
                println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                edge = next_edge
            end
        else
            edge = cw(joint, current_split)
            println("CURRENTLY ON EDGE ($(edge.orig.pos), $(edge.dest.pos))")
            next_face = oppositeface(current_left_face, current_edge)
            if next_face == ll
                nothing
            end
            if current_left_face == hl
                while true
                    # edge.dead = true
                    hideedge(edge)
                    println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                    next_edge = ccwface(current_left_face, edge)
                    if isboundaryedge(next_edge)
                        common_vertex = commonvertex(edge, next_edge)
                        unstickedge(next_edge, common_vertex)
                        if common_vertex == next_edge.dest
                            next_edge.dest = upper_ray.dest
                            next_edge.fr = hl
                            next_edge.ccwd = upper_ray
                            next_edge.cwd = upper_ray.ccwd
                        elseif common_vertex == next_edge.orig
                            next_edge.orig = upper_ray.dest
                            next_edge.fl = hl
                            next_edge.ccwo = upper_ray
                            next_edge.cwo = upper_ray.ccwd
                        end
                        upper_ray.ccwd.ccwd = next_edge
                        upper_ray.cwd = next_edge
                        break
                    end
                    edge = next_edge
                end
            elseif current_left_face == ll
                nothing
            else
                current_endpoints = Set([])
                previous_endpoints = Set([])
                while true
                    current_endpoints = endpoints(edge)
                    inter = intersect(current_endpoints, previous_endpoints)
                    if !isempty(inter)
                        pop!(inter).dead = true
                    end
                    next_edge = ccwface(joint.fr, edge)
                    if edge == next_edge #edge.dead #|| current_left_vertex in endpoints(edge)
                        break
                    end
                    # edge.dead = true
                    hideedge(edge)
                    println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                    edge = next_edge
                    previous_endpoints = current_endpoints
                end
            end
        end

        # joint.fr = current_left_face
        # joint.fl = current_right_face
        # current_left_face.edge = joint
        # current_right_face.edge = joint

        if right_first
            current_right_face = oppositeface(current_right_face, current_edge)
        else
            current_left_face = next_face
        end
        current_vertex = current_split
    end
    return D
end


##
p = [[2,2],[1,1],[3,5]]
q = [[3,3],[5,4],[6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
fixids!(a)
fixids!(b)
checkdcel(a)
checkdcel(b)
test = mergevoronoi(a, b)

##
plotdcel(D, ratio=:equal, sites=true, dead_edges=false, dead_vertices=false)
# checkdcel(test)
