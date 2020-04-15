##
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics

##

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
points = [rand(2) for i in 1:1000]
@time test = voronoihelper(points, io=nothing)

##
@benchmark setdiff(Set([Vertex("1"), Vertex("2")]), Set([Vertex("1"), Vertex("2")]))

##
@benchmark setdiff([Vertex("1"), Vertex("2")], [Vertex("1"), Vertex("2")])


##
@benchmark setdiff((Vertex("1"), Vertex("2")), (Vertex("1"), Vertex("2")))
