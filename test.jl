D = DCEL();
add_join_vertex!(D, [0,0])
add_join_vertex!(D, [1,1], D.vertexlist[1])
add_join_vertex!(D, [1,2], D.vertexlist[2])
add_join_vertex!(D, [3,4], D.vertexlist[2])
add_join_vertex!(D, [1,4], D.vertexlist[3])
add_join_vertex!(D, [0,2], D.vertexlist[1])
joinvertices!(D, D.vertexlist[6], D.vertexlist[3])
joinvertices!(D, D.vertexlist[5], D.vertexlist[4])
joinvertices!(D, D.vertexlist[1], D.vertexlist[3])
joinvertices!(D, D.vertexlist[5], D.vertexlist[6])
add_join_vertex!(D, [2,1], D.vertexlist[2])
joinvertices!(D, D.vertexlist[1], D.vertexlist[7])
addray!(D, D.vertexlist[4], pi/4)
addray!(D, D.vertexlist[5], pi/2)
test = add_join_vertex!(D, [2,4])
splitedge!(D, D.edgelist[7], test)
addray!(D, D.vertexlist[6], 3*pi/4)
test = add_join_vertex!(D, [1,0.5])
splitedge!(D, D.edgelist[11], test)
addray!(D, D.vertexlist[1], -pi/2)
addray!(D, D.vertexlist[7], -pi/2)
addray!(D, D.vertexlist[7], 0)
addray!(D, D.vertexlist[2], pi/8)
split = add_join_vertex!(D, [0,5])
splitedge!(D, D.edgelist[9], split)
addray!(D, D.vertexlist[17], pi/2)
deleteedge!(D, D.edgelist[6])
D.facelist[13].site = [0.5,4]
joinvertices!(D, D.vertexlist[6], D.vertexlist[3])

fixids!(D)
checkdcel(D)
test = plotdcel(D, 1, false)
plot(test)