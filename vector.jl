module Vector
export dot, cross, norm

function cross2d(u::Array, v::Array)::Float64
    return u[1]*v[2] - u[2]*v[1]
end

function dot(u::Array, v::Array)::Float64
    return u[1]*v[1] + u[2]*v[2]
end

function norm(a::Array)::Float64
    return sqrt(sum(a .^2))
end

end
