##
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

##
function mergevoronoi(left::DCEL, right::DCEL, io)
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    fixids!(D)
    push!(D.facelist, fi)
    hl, hr, ll, lr = findextrema(left, right, fi)
    left = right = nothing
    writenothing(io, "TOP FACES: ($hl, $hr)\n")
    writenothing(io, "BOT FACES: ($ll, $lr)\n")

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
    handler.ignore = s
    handler.current_joint = top_ray
    while handler.left_face != ll || handler.right_face != lr
        old_left_face = handler.left_face
        old_right_face = handler.right_face
        foobar(D, handler, io=io)
        if handler.side == :left && old_left_face == hl #by some miracle this also works for the degenerate case
            weldboundary(D, tl, top_ray, :left, io=io)
        elseif handler.side == :right && old_right_face == hr
            weldboundary(D, tr, top_ray, :right, io=io)
        end
    end
    weldboundary(D, bl, bot_ray, :right, io=io)
    updateray!(D, bot_ray, handler.current_vertex, lr, ll)
    weldboundary(D, br, bot_ray, :left, io=io)
    if handler.side == :left
        squeezeedge!(bot_ray, handler.current_vertex, false, previous=handler.current_joint.ccwd, next=handler.current_joint)
    else
        squeezeedge!(bot_ray, handler.current_vertex, false, next=handler.current_joint.cwd, previous=handler.current_joint)
    end

    fi.edge = tl[1]
    hl.edge = tl[1]
    hr.edge = tr[1]
    ll.edge = bl[1]
    lr.edge = br[1]

    filter!(x->!getfield(x, :dead), D.edgelist)
    filter!(x->!getfield(x, :dead), D.vertexlist)
    fixids!(D)
    return D
end


function voronoi(points::Array, io)
    if length(points) == 2
        return voronoitwopoints(points)
    elseif length(points) == 3
        return voronoithreepoints(points)
    else
        split = Int(floor(length(points)/2))
        points_left = points[1:split]
        points_right = points[split+1:end]
        vor_left = voronoi(points_left, io)
        vor_right = voronoi(points_right, io)
        vor = mergevoronoi(vor_left, vor_right, io)
        return vor
    end
    return
end


function voronoihelper(points::Array; io=stdout)
    points = sort(points, by=x->x[1])
    writenothing(io, "POINTS: $points\n")
    out = voronoi(points, io)
    return out
end



##
file = open("log.txt", "w")
# points = [[0.14398304177257648, 0.3135285454600225],
#  [0.14976778967492743, 0.4589056757755534],
#  [0.20297895300489732, 0.5085108793774618],
#  [0.25152768290667016, 0.38248434741059745],
#  [0.2841969985654731, 0.8829357507003006],
#  [0.7369629388643844, 0.6420969230393154]]
# points = temp
# points = [rand(2) for i in 1:6]
points = [[1,1], [2,2], [3,5], [3,3], [5,4], [6,2]]
@time test = voronoihelper(points, io=file)
close(file)
# checkdcel(test, io=nothing)
##
for i in 1:10
    print(i)
    file = open("$i.txt", "w")
    points = [rand(2) for i in 1:1000]
    @time test = voronoihelper(points, io=nothing)
    close(file)
    checkdcel(test, io=nothing)
    plotdcel(test, size=(2000,2000), labels=false, line=:line)
    savefig("$i.png")
end
# plotdcel(test)

##
points = [rand(2) for i in 1:200]
@profiler test = voronoihelper(points, io=nothing)

##
@benchmark intersect(Set([Vertex("1"), Vertex("2")]), Set([Vertex("1"), Vertex("2")]))

##
@benchmark intersect([Vertex("1"), Vertex("2")], [Vertex("1"), Vertex("2")])


##
@benchmark intersect((Vertex("1"), Vertex("2")), (Vertex("1"), Vertex("2")))

##
tangentvector([0,0], 2)
