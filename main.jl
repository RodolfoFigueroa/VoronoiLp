##---
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

##---
function mergevoronoi(left::DCEL, right::DCEL, io)
    fi = mergeinfinitefaces!(left, right)
    global D = joindcel(left, right)
    push!(D.facelist, fi)
    hl, hr, ll, lr = findextrema(left, right, fi)
    left = right = nothing
    # writenothing(io, "TOP FACES: ($hl, $hr)\n")
    # writenothing(io, "BOT FACES: ($ll, $lr)\n")

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
        squeezeedge!(handler.current_vertex, bot_ray, false, previous=handler.current_joint.ccwd, next=handler.current_joint)
    else
        squeezeedge!(handler.current_vertex, bot_ray, false, next=handler.current_joint.cwd, previous=handler.current_joint)
    end

    fi.edge = tl[1]
    hl.edge = tl[1]
    hr.edge = tr[1]
    ll.edge = bl[1]
    lr.edge = br[1]

    filter!(x->!getfield(x, :dead), D.edgelist)
    filter!(x->!getfield(x, :dead), D.vertexlist)
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
    # writenothing(io, "POINTS: $points\n")
    out = voronoi(points, io)
    fixids!(out)
    return out
end



##---
file = open("log.txt", "w")
points = [rand(2) for i in 1:5000]
@profview test = voronoihelper(points, io=file)
close(file)
# fixids!(test)
# checkdcel(test, io=nothing)



##---
# for i in 1:10
#     print(i)
#     file = open("$i.txt", "w")
#     points = [rand(2) for i in 1:1000]
#     @time test = voronoihelper(points, io=nothing)
#     close(file)
#     checkdcel(test, io=nothing)
#     plotdcel(test, size=(2000,2000), labels=false, line=:line)
#     savefig("$i.png")
# end
# # plotdcel(test)

# ##
# t = fullcircle(3.0, 0.01)
# scatter(t[:,1], t[:,2])

# #---
# function testplot(x, y; kwargs...)
#     return plot(x, y; kwargs...)
# end

# testplot([1,2],[3,4], xlims=(1,5))

# #---
# testarray = rand(1000,2)
# testarray2 = [rand(2) for i in 1:1000]
# @time scatter(testarray[:,1], testarray[:,2])
# p = plot()
# @time for i in 1:1000
#     scatter!(p, testarray2[i])
# end

# #---
# p = plotdcel(test, labels=false, line=:solid, sites=true)
