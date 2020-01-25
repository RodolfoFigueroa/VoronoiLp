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

function killface!(edge::Edge, joint::Edge, right::Bool)
    while true
        next_edge = right ? cwface(joint.fl, edge) : ccwface(joint.fr, edge)
        hideedge(edge)
        v = commonvertex(edge, next_edge)
        v.dead = true
        println("KILLED VERTEX: $(v.id)")
        println("KILL1")
        println("KILLED: $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))")
        if edge == next_edge #edge.dead #|| current_left_vertex in endpoints(edge)
            break
        end
        edge = next_edge
    end
    return
end

function firstray!(split::Vertex, angle::Number, hr::Face, hl::Face, right::Bool)
    upper_ray = addray!(D, split, angle+pi, right ? hr : hl, false)
    println("RAY HAS AN ORIGIN $(upper_ray.orig.pos)")
    starter_edge = right ? ccw(upper_ray, split) : cw(upper_ray, split)
    edge = starter_edge
    while true
        hideedge(edge)
        println("KILL2")
        println("KILLED: $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))")
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


function diagramintersection(ir::Array, il::Array, er::Union{Edge,Nothing}, el::Union{Edge,Nothing})
    if ir[2] > il[2]
        println("INTERSECTED RIGHT DIAGRAM FIRST AT $(ir)")
        starter_right = true
        split = createvertex!(D, ir)
        s1, s2 = splitedge!(D, er, split)
    elseif il[2] > ir[2]
        println("INTERSECTED LEFT DIAGRAM FIRST AT $(il)")
        starter_right = false
        split = createvertex!(D, il)
        s1, s2 = splitedge!(D, el, split)
    end
    return starter_right, split, s1, s2
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
    starter_right, split, s1, s2 = diagramintersection(ir, il, er, el)

    current_right_face = current_left_face = Face
    current_right_vertex = current_left_vertex = nothing
    upper_ray = Edge
    if starter_right
        upper_ray = firstray!(split, angle, hr, hl, true)
        current_right_face = oppositeface(hr, upper_ray.cwo)
        current_left_face = hl
        current_right_vertex = split
    else
        upper_ray = firstray!(split, angle, hr, hl, false)
        current_right_face = hr
        current_left_face = oppositeface(hl, upper_ray.ccwo)
        current_left_vertex = split
    end
    current_vertex = split
    weld_set = false
    global weld_vertex = Vertex
    global next_weld_edge = Edge
    finisher_right = Bool
    while true
        angle = perpangle(current_right_face, current_left_face)
        println("\n\n==NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $(current_vertex.pos)==")
        if current_left_face == ll && current_right_face == lr
            # filter!(x -> !getfield(x, :dead), D.edgelist)
            # filter!(x -> !getfield(x, :dead), D.vertexlist)
            lower_ray = addray!(D, current_vertex, angle, finisher_right ? lr : ll, false)
            # return D
            if finisher_right
                lower_ray.ccwd.dead = true
            else
                lower_ray.cwd.dead = true
            end
            if next_weld_edge.orig == weld_vertex
                next_weld_edge.orig = lower_ray.dest
            elseif next_weld_edge.dest == weld_vertex
                next_weld_edge.dest = lower_ray.dest
            else #sanity
                throw("")
            end
            if finisher_right
                squeezeedge!(lower_ray.dest, next_weld_edge, next=lower_ray)
            else
                squeezeedge!(lower_ray.dest, next_weld_edge, previous=lower_ray)
            end
            # fixids!(D)
            return D
        end
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
        current_edge = ir[2]>il[2] ? er : el
        right_first, current_split, s1, s2 = diagramintersection(ir, il, er, el)
        joint = joinvertices!(D, current_vertex, current_split, false)
        joint.fr = current_left_face
        joint.fl = current_right_face
        current_left_face.edge = current_right_face.edge = joint
        if right_first
            edge = joint.ccwd
            println("CURRENTLY ON RIGHT EDGE ($(edge.orig.pos), $(edge.dest.pos))")
            if !weld_set && current_left_face == ll
                println("SPECIAL 1")
                weld_edge = joint.cwo
                while true
                    next_weld_edge = ccwface(lr, weld_edge)
                    if isboundaryedge(next_weld_edge)
                        common_weld_vertex = commonvertex(weld_edge, next_weld_edge)
                        weld_vertex = common_weld_vertex==next_weld_edge.orig ? next_weld_edge.dest : next_weld_edge.orig
                        weld_set = true
                        finisher_right = true
                        break
                    end
                    weld_edge = next_weld_edge
                end
            else
                killface!(edge, joint, true)
            end
        else
            edge = joint.cwd
            println("CURRENTLY ON LEFT EDGE ($(edge.orig.pos), $(edge.dest.pos))")
            if !weld_set && current_right_face == lr
                println("SPECIAL 2")
                weld_edge = joint.ccwo
                while true
                    next_weld_edge = cwface(lr, weld_edge)
                    if isboundaryedge(next_weld_edge)
                        common_weld_vertex = commonvertex(weld_edge, next_weld_edge)
                        weld_vertex = common_weld_vertex==next_weld_edge.orig ? next_weld_edge.dest : next_weld_edge.orig
                        weld_set = true
                        finisher_right = false
                        break
                    end
                    weld_edge = next_weld_edge
                end
                killface!(edge, joint, false)
            elseif current_left_face == hl
                while true
                    # edge.dead = true
                    println("KILL3")
                    println("KILLED: $(edge.id) ($(edge.orig.pos), $(edge.dest.pos))")
                    hideedge(edge)
                    next_edge = ccwface(current_left_face, edge)
                    commonvertex(edge, next_edge).dead = true
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
            else
                killface!(edge, joint, false)
            end
        end

        # joint.fr = current_left_face
        # joint.fl = current_right_face
        # current_left_face.edge = joint
        # current_right_face.edge = joint

        if right_first
            current_right_face = oppositeface(current_right_face, current_edge)
        else
            current_left_face = oppositeface(current_left_face, current_edge)
        end
        current_vertex = current_split
        # counter +=1
        # if counter==2
        #     return D
        # end
    end
    return D
end


##
p = [[2,2],[1,1],[3,5]]
q = [[3,3],[5,4],[6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
fixids!(a)
# fixids!(b)
checkdcel(a)
checkdcel(b)
test = mergevoronoi(a, b)

##
plotdcel(test, ratio=:equal, sites=true, dead_edges=true, dead_vertices=false)
# checkdcel(test)
