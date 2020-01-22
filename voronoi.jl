module DVoronoi
using Statistics

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