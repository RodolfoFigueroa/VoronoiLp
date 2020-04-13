##
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

##
function voronoi(start_points::Array, io=stdout)
    points = sort(start_points, by=x->x[1])
    if length(points) == 2
        return voronoitwopoints(points)
    elseif length(points) == 3
        return voronoithreepoints(points)
    else
        split = Int(floor(length(points)/2))
        points_left = points[1:split]
        points_right = points[split+1:end]
        vor_left = voronoi(points_left)
        vor_right = voronoi(points_right)
        # return points_left, points_right

        # return vor_left, vor_right
        vor = mergevoronoi(vor_left, vor_right, io)
        return vor
    end
    return
end

function mergevoronoi(left::DCEL, right::DCEL, io)
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    fixids!(D)
    push!(D.facelist, fi)
    hl, hr, ll, lr = findextrema(left, right)
    write(io, "TOP FACES: ($hl, $hr)\n")
    write(io, "BOT FACES: ($ll, $lr)\n")

    tl, tr, bl, br = openface.([hl, hr, ll, lr], [:ccw, :cw, :cw, :ccw], io=io)
    starter_left = getfield(tl[1], tl[2])
    starter_right = getfield(tr[1], tr[2])
    top_ray = createfloatingedge!(D)
    bot_ray = createfloatingedge!(D) #comment out if you cannot print the edgelist

    handler = Handler(hl, hr, starter_left, starter_right)
    right_first, split, s = highestintersection(D, handler, io=io)
    if right_first == :right
        edge = cwface(hr, s[1]) == s[2] ? s[2] : s[1]
        t = tr
    else
        edge = ccwface(hl, s[1]) == s[2] ? s[2] : s[1]
        t = tl
    end
    weldboundary(D, t, top_ray, right_first, io=io)
    updateray!(D, top_ray, split, hl, hr)
    updatehandler!(handler, edge, split, right_first, starter_left, starter_right, io=io)
    handler.ignore = [] #s
    handler.current_joint = top_ray
    bottom_vertex = nothing
    counter = 0
    while handler.left_face != ll || handler.right_face != lr
        write(io, "\n%%%TOP METHOD%%%\n")
        old_left_face = handler.left_face
        old_right_face = handler.right_face
        foobar(D, handler, false, io=io)
        old_side = handler.side
        println(old_side)
        println(old_left_face)
        println(old_right_face)
        if old_side == :left
            if old_left_face == hl
                println("T1")
                weldboundary(D, tl, top_ray, :left, io=io)
            end
        else
            if old_right_face == hr
                println("T3")
                weldboundary(D, tr, top_ray, :right, io=io)
            end
        end

        counter += 1
        if counter == 1
            # return
        end
    end
    weldboundary(D, bl, bot_ray, :right)
    updateray!(D, bot_ray, handler.current_vertex, lr, ll)
    weldboundary(D, br, bot_ray, :left)
    if handler.side == :left
        squeezeedge!(handler.current_vertex, bot_ray, false, previous=handler.current_joint.ccwd, next=handler.current_joint)
    else
        # settopology!(bot_ray, cwo=handler.current_joint, ccwo=handler.current_joint.cwd)
        squeezeedge!(handler.current_vertex, bot_ray, false, next=handler.current_joint.cwd, previous=handler.current_joint)
    end
    fi.edge = bot_ray.ccwd

    hl.edge = tl[1]
    hr.edge = tr[1]
    ll.edge = bl[1]
    lr.edge = br[1]

    filter!(x->!getfield(x, :dead), D.edgelist)
    filter!(x->!getfield(x, :dead), D.vertexlist)
    fixids!(D)
    return D
end
##
file = open("log.txt", "w")
# points = temp
# points = [[0.14398304177257648, 0.3135285454600225],
#  [0.14976778967492743, 0.4589056757755534],
#  [0.20297895300489732, 0.5085108793774618],
#  [0.25152768290667016, 0.38248434741059745],
#  [0.2841969985654731, 0.8829357507003006],
#  [0.7369629388643844, 0.6420969230393154]]
points = [[1,1], [2,2], [3,5], [3,3], [5,4], [6,2]]
# points = [rand(2) for i in 1:6]
test = voronoi(points, file)
close(file)
checkdcel(D)

##
for i in 1:500
    global file = open("log.txt", "w")
    global points = [rand(2) for i in 1:6]
    test = voronoi(points, file)
    close(file)
    checkdcel(test)
end
# plotdcel(test)
