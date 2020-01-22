##
include("./vector.jl")
include("./struct.jl")
using .DVector
using .DStruct
using Plots
using Statistics
##
# D = DCEL();
# add_join_vertex!(D, [0,0])
# add_join_vertex!(D, [1,1], D.vertexlist[1])
# add_join_vertex!(D, [1,2], D.vertexlist[2])
# add_join_vertex!(D, [3,4], D.vertexlist[2])
# add_join_vertex!(D, [1,4], D.vertexlist[3])
# add_join_vertex!(D, [0,2], D.vertexlist[1])
# joinvertices!(D, D.vertexlist[6], D.vertexlist[3])
# joinvertices!(D, D.vertexlist[5], D.vertexlist[4])
# joinvertices!(D, D.vertexlist[1], D.vertexlist[3])
# joinvertices!(D, D.vertexlist[5], D.vertexlist[6])
# add_join_vertex!(D, [2,1], D.vertexlist[2])
# joinvertices!(D, D.vertexlist[1], D.vertexlist[7])
# addray!(D, D.vertexlist[4], pi/4)
# addray!(D, D.vertexlist[5], pi/2)
# test = add_join_vertex!(D, [2,4])
# splitedge!(D, D.edgelist[7], test)
# addray!(D, D.vertexlist[6], 3*pi/4)
# test = add_join_vertex!(D, [1,0.5])
# splitedge!(D, D.edgelist[11], test)
# addray!(D, D.vertexlist[1], -pi/2)
# addray!(D, D.vertexlist[7], -pi/2)
# addray!(D, D.vertexlist[7], 0)
# addray!(D, D.vertexlist[2], pi/8)
# split = add_join_vertex!(D, [0,5])
# splitedge!(D, D.edgelist[9], split)
# addray!(D, D.vertexlist[17], pi/2)
# deleteedge!(D, D.edgelist[6])
# D.facelist[13].site = [0.5,4]
# joinvertices!(D, D.vertexlist[6], D.vertexlist[3])
#
# fixids!(D)
# checkdcel(D)
# test = plotdcel(D, 1, false)
# plot(test)
##

leftofedge(e::Edge, v::Vertex)::Bool = leftofedge(e, v.pos)

function perpangle(start::Array, finish::Array)::Float64
    return atan(finish[2]-start[2], finish[1]-start[1]) + pi/2
end

function perpangle(f1::Face, f2::Face)::Float64
    return perpangle(f1.site, f2.site)
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
        return u>=0 && t>=0 ? p+t*r : [NaN,NaN]
    else
        return u>=0 && 0<=t<=1 ? p+t*r : [NaN,NaN]
    end
    return
end

function edgerayintersection(edge::Edge, start::Array, finish::Array)::Array
    p = edge.orig.pos
    r = (edge.dest.pos-p)/norm(edge.dest.pos-p)
    s = [cosatan(finish[2]-start[2],finish[1]-start[1]),sinatan(finish[2]-start[2],finish[1]-start[1])]
    num = start - p
    den = cross2d(r, s)
    t = cross2d(num, s)/den
    u = cross2d(num, r)/den
    return u>=0 && 0<=t<=1 ? p+t*r : [NaN,NaN]
end

function facerayintersection(f::Face, start::Array, angle::Number, starter_edge::Union{Edge,Nothing}=nothing)
    edge = isnothing(starter_edge) ? f.edge : starter_edge
    while true
        if !isboundaryedge(edge)
            inter = edgerayintersection(edge, start, angle)
            if (!isnan(inter[1])) && (isnothing(starter_edge) || edge != starter_edge)
                return edge, inter
            end
        end
        edge = ccwface(f, edge)
    end
    return
end

function oppositeface(f::Face, e::Edge)::Face
    println("CHECKING OPPOSITE FACE OF: ")
    show(f)
    println("WITH MIRROR EDGE: $(e.orig.pos), $(e.dest.pos)")
    println("MIRROR FACE RIGHT: $(e.fr)")
    println("MIRROR FACE LEFT: $(e.fl)")
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

function mergevoronoi(left::DCEL, right::DCEL)
    D = joindcel(left, right)
    hl, ll = findextrema(left)
    hr, lr = findextrema(right)
    angle = perpangle(hr, hl)
    println("INITIAL RAY HAS AN ANGLE: $(rad2deg(angle))")
    el, il = facerayintersection(hl, midpoint(hr, hl), angle)
    er, ir = facerayintersection(hr, midpoint(hr, hl), angle)
    split = Vertex
    starter_right = Bool
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
    current_right = Face
    current_left = Face
    ray = Edge
    if starter_right
        ray = addray!(D, split, angle+pi, hr)
        println("RAY HAS AN ORIGIN $(ray.orig.pos)")
        println("RAY CWO: $(ray.cwo.orig.pos), $(ray.cwo.dest.pos)")
        println("RAY CCWO: $(ray.ccwo.orig.pos), $(ray.ccwo.dest.pos)")
        println("RAY FACE RIGHT: $(ray.fr)")
        println("RAY FACE LEFT: $(ray.fl)")
        # starter_edge = cw(s1,split) == ray ? s1 : s2
        starter_edge = ccw(ray, split)
        edge = starter_edge
        while true
            next = cwface(ray.fl, edge)
            # deleteedge!(D, edge)
            edge.dead = true
            if isboundaryedge(edge)
                break
            end
            commonvertex(edge, next).dead = true
            edge = next
        end
        current_right = oppositeface(ray.fl, starter_edge)
        current_left = hl
        current_vertex = split
        ray.fl = hl
        ray.fr = hr
    end
    while true
        if current_left == ll && current_right == lr
            return
        end
        angle = perpangle(current_right, current_left)
        println("NEXT BISECTOR HAS AN ANGLE: $(rad2deg(angle))")
        el, il = facerayintersection(current_left, current_vertex.pos, angle)
        er, ir = facerayintersection(current_left, current_vertex.pos, angle)
        current_edge = Edge
        current_split = Vertex
        right_first = Bool
        if ir[2] > il[2]
            println("BISECTOR INTERSECTED RIGHT EDGE: $(er.orig.pos) $(er.dest.pos)")
            println("AT: $(ir)")
            right_first = true
            current_split = createvertex!(D, ir)
            current_edge = er
        else
            println("BISECTOR INTERSECTED LEFT EDGE $(el.orig.pos) $(el.dest.pos)")
            println("AT: $(il)")
            right_first = false
            current_split = createvertex!(D, il)
            current_edge = el
        end
        s1, s2 = splitedge!(D, current_edge, current_split)
        joint = joinvertices!(D, current_vertex, current_split)
        if right_first
            edge = cw(s1,split) == joint ? s1 : s2
            while true
                edge.dead = true
                if split in endpoints(edge)
                    break
                end
                edge = cwface(joint.fl, edge)
            end
        else
            edge = ccw(s1,split) == joint ? s1 : s2
            while true
                edge.dead = true
                if split in endpoints(edge)
                    break
                end
                edge = ccwface(joint.fr, edge)
            end
        end

        joint.fr = current_left
        joint.fl = current_right
        current_left.edge = joint
        current_right.edge = joint

        if right_first
            current_right = oppositeface(current_right, current_edge)
        else
            current_left = oppositeface(current_left, current_edge)
        end

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
# checkdcel(test)
