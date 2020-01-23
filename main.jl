##
include("./vector.jl")
include("./struct.jl")
using .DVector
using .DStruct
using Plots
using Statistics

##
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
    a, b, c = pointccw(points) ? bisectorangles(v, u, w) : bisectorangles(u, v, w)
    D = DCEL()
    m = createvertex!(D, circlethreepoints(u, v, w))
    p = addray!(D, m, a)
    r = addray!(D, m, c)
    # return D
    outeredge = angleccw(a,b,c) ? p.cwd : r.cwd
    dummy = createdummyvertex!(D, m, b)
    s1, s2 = splitedge!(D, outeredge, dummy)
    joinvertices!(D, m, dummy, p2=s2)
    p.fl.site = v
    p.fr.site = w
    p.ccwo.fl.site = u
    return D
end

function findextrema(D::DCEL)
    positions = getfield.(D.facelist, :site)
    out = [-Inf, -Inf]
    index_high = 1
    index_low = 1
    for i in 1:length(positions)
        if !isnothing(positions[i]) && positions[i][2] > positions[index_high][2]
            index_high = i
        end
        if !isnothing(positions[i]) && positions[i][2] < positions[index_low][2]
            index_low = i
        end
    end
    return D.facelist[index_high], D.facelist[index_low]
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

function edgerayintersection(edge::Edge, q::Array, angle::Float64)::Array
    p = edge.orig.pos
    r = (edge.dest.pos-p)/norm(edge.dest.pos-p)
    s = [cos(angle),sin(angle)]
    num = q - p
    den = cross2d(r, s)
    t = cross2d(num, s)/den
    u = cross2d(num, r)/den
    if isstrutedge(edge)
        return u>=0 && 0<=t #=&& cross2d(p-q, r) != 0=# ? p+t*r : [NaN,NaN]
    else
        return u>=0 && 0<=t<=1 ? p+t*r #=&& cross2d(p-q, r) != 0=# : [NaN,NaN]
    end
    return
end

# function edgerayintersection(edge::Edge, start::Array, finish::Array)::Array
#     p = edge.orig.pos
#     r = (edge.dest.pos-p)/norm(edge.dest.pos-p)
#     s = [cosatan(finish[2]-start[2],finish[1]-start[1]),sinatan(finish[2]-start[2],finish[1]-start[1])]
#     num = start - p
#     den = cross2d(r, s)
#     t = cross2d(num, s)/den
#     u = cross2d(num, r)/den
#     return u>=0 && 0<=t<=1 ? p+t*r : [NaN,NaN]
# end

function facerayintersection(f::Face, start::Array, angle::Number, clockwise::Bool, ignore::Array=[])
    starter_edge = f.edge
    edge = starter_edge
    println("\nCHECKING INTERSECTION OF FACE $(f)")
    while true
        # sleep(0.1)
        println("\nCHECKING EDGE ($(edge.orig.pos), $(edge.dest.pos))")
        println("WITH LEFT FACE: $(edge.fl) AND RIGHT FACE: $(edge.fr)")
        if !isboundaryedge(edge) && !(edge in ignore)
            inter = edgerayintersection(edge, start, angle)
            if !isnan(inter[1])
                return edge, inter
            end
        end
        edge = clockwise ? cwface(f,edge) : ccwface(f,edge)
        if edge == starter_edge
            return nothing, [-Inf,-Inf]
        end
    end
    return
end

function oppositeface(f::Face, e::Edge)::Face
    println("CHECKING OPPOSITE FACE OF: $(f)")
    println("\nWITH MIRROR EDGE: $(e.orig.pos), $(e.dest.pos)")
    println("MIRROR RIGHT FACE: $(e.fr)")
    println("MIRROR LEFT FACE: $(e.fl)")
    if e.fr.site == f.site
        println("FOUND OPPOSITE: ")
        show(e.fl)
        return e.fl
    elseif e.fl.site == f.site
        println("FOUND OPPOSITE: ")
        show(e.fr)
        return e.fr
    else
        throw("Face $(f) has a topology problem")
    end
    return
end

#When I wrote this, only God and I understood what was going on.
#Now, only God knows.
function mergevoronoi(left::DCEL, right::DCEL)
    fi = mergeinfinitefaces!(left, right)
    D = joindcel(left, right)
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
        upper_ray = addray!(D, split, angle+pi, hr)
        println("RAY HAS AN ORIGIN $(ray.orig.pos)")
        starter_edge = ccw(upper_ray, split)
        edge = starter_edge
        while true
            next = cwface(upper_ray.fl, edge)
            edge.dead = true
            println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
            if isboundaryedge(edge)
                break
            end
            commonvertex(edge, next).dead = true
            edge = next
        end
        current_right_face = oppositeface(upper_ray.fl, starter_edge)
        current_left_face = hl
        current_vertex = split
        upper_ray.fl = hl
        upper_ray.fr = hr
        current_right_vertex = current_vertex
        welded = false
    end
    while true
        if current_left_face == ll && current_right_face == lr
            return D
        end
        angle = perpangle(current_right_face, current_left_face)
        println("\n\n==NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle)) AND STARTS AT: $(current_vertex.pos)==")
        er, ir = facerayintersection(current_right_face, current_vertex.pos, angle, true, [s1,s2])
        if !isnothing(er)
            println("\n=BISECTOR INTERSECTED RIGHT EDGE: ($(er.orig.pos), $(er.dest.pos) AT: $(ir)=")
        end
        el, il = facerayintersection(current_left_face, current_vertex.pos, angle, false, [s1,s2])
        if !isnothing(el)
            println("\n=BISECTOR INTERSECTED LEFT EDGE ($(el.orig.pos), $(el.dest.pos) AT: $(il)=")
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
        joint = createedge!(D, current_vertex, current_split)
        squeezeedge!(current_vertex, joint)
        squeezeedge!(current_split, joint)
        joint.fr = current_left_face
        joint.fl = current_right_face
        current_left_face.edge = joint
        current_right_face.edge = joint

        if right_first
            edge = ccw(joint, current_split)
            while true
                if edge.dead || current_right_vertex in endpoints(edge)
                    break
                end
                edge.dead = true
                println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                edge = cwface(joint.fl, edge)
            end
        else
            edge = cw(joint, current_split)
            println("CURRENTLY ON EDGE ($(edge.orig.pos), $(edge.dest.pos))")
            if !welded
                @assert current_left_face == hl
                while true
                    if isboundaryedge(edge)
                        if edge.orig.pos[1] < edge.dest.pos[1]
                            println("THE WELD VERTEX IS TO THE LEFT OF THE STARTING RAY")
                            show(edge)
                            edge.dest = upper_ray.dest
                            edge.fr = hl
                            edge.ccwd = upper_ray
                            edge.cwd = upper_ray.ccwd
                        else
                            println("THE WELD VERTEX IS TO THE RIGHT OF THE STARTING RAY")
                            show(edge)
                            edge.orig = upper_ray.dest
                            edge.fl = hl
                            edge.ccwo = upper_ray
                            edge.cwo = upper_ray.ccwd
                        end
                        welded = true
                        break
                    end
                    edge.dead = true
                    println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                    edge = ccwface(current_left_face, edge)
                end
            else
                while true
                    if edge.dead || current_left_vertex in endpoints(edge)
                        break
                    end
                    edge.dead = true
                    println("KILLED: ($(edge.orig.pos), $(edge.dest.pos))")
                    edge = ccwface(joint.fr, edge)
                end
            end
        end

        joint.fr = current_left_face
        joint.fl = current_right_face
        current_left_face.edge = joint
        current_right_face.edge = joint

        if right_first
            current_right_face = oppositeface(current_right_face, current_edge)
        else
            current_left_face = oppositeface(current_left_face, current_edge)
        end
        current_vertex = current_split

    end
    return D
end

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

function mergeinfinitefaces!(a::DCEL, b::DCEL)
    fi = Face("i")
    resetfacelist!(a, fi)
    resetfacelist!(b, fi)
    return fi
end

function joindcel(a::DCEL, b::DCEL)::DCEL
    vertexlist = vcat(a.vertexlist, b.vertexlist)
    edgelist = vcat(a.edgelist, b.edgelist)
    facelist = vcat(a.facelist, b.facelist)
    dummylist = vcat(a.dummylist, b.dummylist)
    return DCEL(edgelist, vertexlist, facelist, dummylist)
end

function circletwopointsradius(p::Array, q::Array, r::Number)
    x1, y1 = p
    x2, y2 = q
    q = sqrt((x2-x1)^2 + (y2-y1)^2)
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
plotdcel(test, ratio=:equal, sites=true)
# checkdcel(test)
