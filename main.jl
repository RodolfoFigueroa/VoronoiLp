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
        println("LEFT POINTS: $points_left")
        println("RIGHT POINTS: $points_right")
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
    b1 = hl == ll
    b2 = hr == lr
    starter_left = getfield(tl[1], tl[2])
    starter_right = getfield(tr[1], tr[2])

    helper1(tl)
    helper1(tr)
    helper1(bl)
    helper1(br)

    top_dummy = createvertex!(D, [0,0])
    bot_dummy = createvertex!(D, [0,0])
    top_dummy.original = bot_dummy.original = false
    println("READY")
    println("TOP LEFT VERTEX: $(getfield(tl[1],tl[3]))")
    println("TOP RIGHT VERTEX: $(getfield(tr[1],tr[3]))")
    println("BOT LEFT VERTEX: $(getfield(bl[1],bl[3]))")
    println("BOT RIGHT VERTEX: $(getfield(br[1],br[3]))")
    println("TOP LEFT EDGE: $(tl[1])")
    println("TOP RIGHT EDGE: $(tr[1])")
    println("BOT LEFT EDGE: $(bl[1])")
    println("BOT RIGHT EDGE: $(br[1])")
    top_left_joint = joinvertices!(D, top_dummy, getfield(tl[1],tl[3]), false)
    top_right_joint = joinvertices!(D, top_dummy, getfield(tr[1],tr[3]), false)
        bot_left_joint =
    
    handler = Handler(hl, hr, starter_left, starter_right)
    right_first, split, s = highestintersection(D, handler, io=io)
    if right_first == :right
        edge = cwface(hr, s[1]) == s[2] ? s[2] : s[1]
    else
        edge = ccwface(hl, s[1]) == s[2] ? s[2] : s[1]
    end
    updatehandler!(handler, edge, split, right_first, starter_left, starter_right, io=io)
    handler.ignore = s

    top_vertex = handler.current_vertex
    top = createray!(D, top_vertex, hl, hr)
    top_ray = top.edge
    hl.edge = hr.edge = top_ray
    handler.current_joint = top_ray

    return D

    if !b1
        # helper1(tl)
        j = helper2(D, top, tl, top_ray, :left)
        j.fl = hl
        j.fr = fi
        fi.edge = j
    else
        tl[1].dead = true
        bl[1].dead = true
        println("{{{{}}}}")
        println(getfield(bl[1], bl[3]))
        joinvertices!(D, top, getfield(bl[1], bl[3]), false)
    end
    if !b2
        # helper1(tr)
        j = helper2(D, top, tr, top_ray, :right)
        j.fl = hr
        j.fl = fi
    end
    bottom_vertex = nothing
    counter = 0
    while true
        if handler.current_left_face == ll && handler.current_right_face == lr
            bottom_vertex = handler.current_vertex
            break
        end
        write(io, "\n%%%TOP METHOD%%%\n")
        foobar(D, handler, io=io)
        counter += 1
        if counter == 2
            return
        end
    end
    return D
    bot = createray!(D, bottom_vertex, lr, ll)
    if handler.side == :right
        bot.edge.cwo = handler.current_joint
        handler.current_joint.ccwd = bot.edge
    else
        bot.edge.ccwo = handler.current_joint
        handler.current_joint.cwd = bot.edge
    end
    bot_ray = bot.edge
    ll.edge = lr.edge = bot_ray
    if !b1
        # helper1(bl)
        j = helper2(D, bot, bl, bot_ray, :right)
        j.fr = ll
        j.fl = fi
    else
        tl[1].dead = true
        bl[1].dead = true
        joint = joinvertices!(D, top, bot, false, no=top_ray, pd=bot_ray)
        joint.fl = hl
        joint.fr = fi
        fi.edge = joint
    end
    if !b2
        # helper1(br)
        j = helper2(D, bot, br, bot_ray, :left)
        j.fl = lr
        j.fr = fi
    else
        tr[1].dead = true
        br[1].dead = true
        joint = joinvertices!(D, top, bot, false, po=top_ray, nd=bot_ray)
        joint.fl = fi
        joint.fr = hr
        fi.edge = joint
    end
    filter!(x->!getfield(x, :dead), D.edgelist)
    filter!(x->!getfield(x, :dead), D.vertexlist)
    fixids!(D)
    return D
end
##
file = open("log.txt", "w")
# points = temp
points = [[0.14398304177257648, 0.3135285454600225],
 [0.14976778967492743, 0.4589056757755534],
 [0.20297895300489732, 0.5085108793774618],
 [0.25152768290667016, 0.38248434741059745],
 [0.2841969985654731, 0.8829357507003006],
 [0.7369629388643844, 0.6420969230393154]]
test = voronoi(points, file)
close(file)

##
for i in 1:500
    global file = open("log.txt", "w")
    global points = [rand(2) for i in 1:6]
    test = voronoi(points, file)
    close(file)
    checkdcel(test)
end
# plotdcel(test)
