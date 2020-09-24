##---
using Revise
includet("./struct.jl")
using .DStruct
using Plots
using Statistics
using Profile, PProf, BenchmarkTools

##---
file = open("log.txt", "w")
points = [rand(2) for i in 1:10000]
dumpoints = [[1,1], [3,5], [0,2]]

##---
file = open("log.txt", "w")
@time test = voronoihelper(points, io=file)
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