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
    if !b1
        j = helper(D, top, top_ray, tl, :left)
        j.fr = fi
        fi.edge = j
    end
    if !b2
        j = helper(D, top, top_ray, tr, :right)
        j.fl = fi
    end
    handler.current_joint = top_ray
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
        if counter == 1
            # return
        end
    end
    bot = createray!(D, bottom_vertex, lr, ll)
    if handler.side == :right
        bot.edge.cwo = handler.current_joint
        handler.current_joint.ccwd = bot.edge
    else
        bot.edge.ccwo = handler.current_joint
        handler.current_joint.cwd = bot.edge
    end
    return D
    bot_ray = bot.edge
    if !b1
        j = helper(D, bot, bot_ray, bl, :right)
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
        j = helper(D, bot, bot_ray, br, :left)
        j.fr = fi
    else
        tr[1].dead = true
        br[1].dead = true
        joint = joinvertices!(D, top, bot, false, po=top_ray, nd=bot_ray)
        joint.fl = fi
        joint.fr = hr
        fi.edge = joint
    end
    hl.edge = hr.edge = top_ray
    ll.edge = lr.edge = bot_ray
    filter!(x->!getfield(x, :dead), D.edgelist)
    filter!(x->!getfield(x, :dead), D.vertexlist)
    fixids!(D)
    return D
end
##
file = open("log.txt", "w")
# p = [[0.5708646705995277, 0.8481649641678017],
#     [0.6288045968631606, 0.24546568227328835],
#     [0.6500802866529276, 0.21643130460624627]]
# q = [[0.9464809834988424, 0.3435419772758477],
#     [0.7576321201041849, 0.611379875342847],
#  [0.7310684327092163, 0.015002978227499186]]
p = [[1,1], [2,2], [3,5]]
q = [[3,3], [5,4], [6,2]]
a = voronoithreepoints(p)
b = voronoithreepoints(q)
checkdcel(a)
checkdcel(b)
##
test = mergevoronoi(a,b, file)
close(file)


##
file = open("log.txt", "w")
points = [rand(2) for i in 1:6]
test = voronoi(points, file)
close(file)
plotdcel(test)


##
plotdcel(test, scale=5)

##
checkdcel(test)

##



##
t = voronoi(vcat(q,p))
