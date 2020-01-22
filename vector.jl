module DVector
export dot, cross2d, norm, distance, pointccw, midpoints, angleccw

function cross2d(u::Array, v::Array)::Float64
    return u[1]*v[2] - u[2]*v[1]
end

function dot(u::Array, v::Array)::Float64
    return u[1]*v[1] + u[2]*v[2]
end

function norm(a::Array)::Float64
    return sqrt(sum(a .^2))
end

function distance(p::Array, q::Array)::Float64 
	return norm(p .-q)
end

function cosatan(y::Number, x::Number)::Float64
    return x/hypot(x,y)
end

function sinatan(y::Number, x::Number)::Float64
    if x==0
        if y>0
            return 1
        elseif y<0
            return -1
        else
            return 0
        end
    else
        return y/hypot(x,y)
    end
    return
end

function pointccw(array::Array)::Bool
    sum = 0
    for i in 2:length(array)
        sum += (array[i][1] - array[i-1][1])*(array[i][2] + array[i-1][2])
    end
    sum += (array[1][1] - array[end][1])*(array[1][2] + array[end][2])
    return sum <= 0
end

function midpoints(u::Array, v::Array, w::Array)::Tuple 
    return mean([u,v]), mean([v,w]), mean([w,u])
end

function angleccw(a::Number, b::Number, c::Number)::Bool
    return sin(a-b) - sin(a-c) + sin(b-c) <= 0
end

end
