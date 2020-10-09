##---
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics
using Profile, PProf, BenchmarkTools

##---
file = open("log.txt", "w")
points = [rand(2) for i in 1:6]
dumpoints = [[1,1], [3,5], [0,2]]

##---
file = open("log.txt", "w")
voronoihelper(dumpoints, io=file)
@time voronoihelper(points, io=file)
# fixids!(test)
# checkdcel(test, io=nothing)
close(file)
# checkdcel(test)

##---
reset()
p1 = [[-1,1], [-2,2], [-3,4]]
p2 = [[1,1], [3,5]]
p = vcat(p1, p2)

file = open("log.txt", "w")
test = voronoihelper(p; io=file)
close(file)

##---
@time commonvertex(v1, v2, v3)
@time commonvertex2(v1, v2, v3) 