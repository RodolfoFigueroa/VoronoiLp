module DVector
export dot, cross2d, norm, distance, pointccw, midpoints, angleccw, cosatan, sinatan, circlethreepoints

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

function circletwopointsradius(p::Array, q::Array, r::Number)
    x1, y1 = p
    x2, y2 = q
    q = sqrt((x2-x1)^2 + (y2-y1)^2)
    y3 = (y1+y2)/2
    x3 = (x1+x2)/2
    basex = sqrt(r^2 - (q/2)^2) * (y1-y2)/q
    basey = sqrt(r^2 - (q/2)^2) * (x2-x1)/q
    return [x3+basex, y3+basey], [x3-basex, y3-basey]
end

function circlethreepoints(p::Array, q::Array, r::Array)
    x1,y1 = p
    x2,y2 = q
    x3,y3 = r
    den = 2*(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))
    x = (x1^2+y1^2)*(y3-y2)+(x2^2+y2^2)*(y1-y3)+(x3^2+y3^2)*(y2-y1)
    y = (x1^2+y1^2)*(x2-x3)+(x2^2+y2^2)*(x3-x1)+(x3^2+y3^2)*(x1-x2)
    return [x/den,y/den]
end

end
