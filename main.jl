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
voronoihelper(dumpoints, io=file)
test = voronoihelper(points, io=file)
close(file)
# fixids!(test)
# checkdcel(test, io=nothing)