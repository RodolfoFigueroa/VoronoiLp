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
Profile.clear()
@time test = voronoihelper(points, io=file)
close(file)
# fixids!(test)
# checkdcel(test, io=nothing)

##---
@benchmark voronoihelper(points) setup=(points=[rand(2) for i in 1:5000]) seconds=20


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

t = rand(100000)
function a1(x)
    return findfirst(i->i>0.5, x)
end

function a2(x)
    for i in x
        if i >0.5
            return i
        end
    end
    return
end

@time a1(t)
@time a2(t)