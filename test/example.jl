# code from the example notebook for automatic testing

using PyPlot

using Coils
using Coils.CoilsPlot

g, vertex_positions = cuboid_system([1, 1, 1], [3, 3, 3])

figure(figsize = (6, 6))
plot_vertices(vertex_positions)
plot_edges(g, vertex_positions, standalone = false)

poi = cuboid_poi([0.75, 0.75, 0.75], [0.0, 0.0, 0.0], [10, 10, 10], filled = false)

figure(figsize = (6, 6))
plot_system(g, vertex_positions, poi)

Bgoal(x) = [0, 100e-6, 0]

cells = find_cells(g)

figure(figsize = (8, 8))
plot_vertices(vertex_positions, labels = false)
plot_cells(cells, vertex_positions)

M = system_matrix(poi, vertex_positions, cells)

Bpoi = vcat(Bgoal.(poi)...)

optI = M \ Bpoi

figure(figsize = (10, 10))
plot_vertices(vertex_positions, labels = false)
current_norm = 1000 / maximum(abs, optI)
plot_cells(cells, vertex_positions, optI * current_norm)

edgecurrents = find_edgecurrents(g, cells, optI)

figure(figsize = (8, 8))
plot_vertices(vertex_positions, labels = false)
plot_edge_currents(g, edgecurrents, vertex_positions)

simpleloops, simpleloopscurrents = find_all_simpleloops(g, edgecurrents)

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[1], vertex_positions, color = "orange")
plot_loop(simpleloops[2], vertex_positions, color = "orange")
title("Current: $(simpleloopscurrents[1] * current_norm)")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[3], vertex_positions, color = "green")
plot_loop(simpleloops[4], vertex_positions, color = "green")
title("Current: $(simpleloopscurrents[3] * current_norm)")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[5], vertex_positions, color = "blue")
plot_loop(simpleloops[6], vertex_positions, color = "blue")
title("Current: $(simpleloopscurrents[5] * current_norm)")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[7], vertex_positions, color = "indigo")
plot_loop(simpleloops[8], vertex_positions, color = "indigo")
title("Current: $(simpleloopscurrents[7] * current_norm)")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
plot_loop(simpleloops[9], vertex_positions, color = "grey")
plot_loop(simpleloops[10], vertex_positions, color = "grey")
title("Current: $(simpleloopscurrents[9] * current_norm)")

figure()
plot_deviation_histogram(poi, simpleloops, simpleloopscurrents, vertex_positions, Bgoal)

figure()
plot_loops_field(cells, optI, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-0.5, 0.5], Bref = Bgoal, levels = -10:10)

g, vertex_positions = cuboid_system([1, 2, 1], [3, 6, 3],
    skipfaces = [false, false, true, true, false, false])
poi = cuboid_poi([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [10, 10, 10], filled = true)

# the wires used to construct the coil
elemcurrents = [10.0, 1.0, 0.1]

figure(figsize = (6, 6))
plot_system(g, vertex_positions, poi)

# Find out the location of the openings
extrema(getindex.(vertex_positions, 2))

# prepare the large cells 
front_face_unordered = [ i for i in eachindex(vertex_positions) if vertex_positions[i][2] < -0.9 ]
front_face = order_cell(g, front_face_unordered)
            
back_face_unordered = [ i for i in eachindex(vertex_positions) if vertex_positions[i][2] > 0.9 ]
back_face = order_cell(g, back_face_unordered)
                        
initialcells = [front_face, back_face]

Bgoal(x) = [√2 * 100e-6 * x[2], √2 * 100e-6 * x[1], 0]

simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal,
    verbose = true, initialcells = [front_face, back_face],
    mincurrent = minimum(elemcurrents) / 2)

simpleloopscurrents_decomp = decompose_currents(simpleloopscurrents, elemcurrents)
simpleloopscurrents_real = real_decomposed_currents(simpleloopscurrents_decomp, elemcurrents)

figure()
plot_deviation_histogram(poi, simpleloops, simpleloopscurrents_real, vertex_positions, Bgoal)

figure()
plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-1, 1], Bref = Bgoal, levels = -20:4:20)

figure()
plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-1, 1], levels = -50:5:50)

for i in eachindex(simpleloopscurrents)
    print("#$(lpad(i,2)): ")
    print("$(lpad(simpleloopscurrents[i], 20)) = ")
    print("$(lpad(simpleloopscurrents_decomp[i][1], 2)) * 10A  +  ")
    print("$(simpleloopscurrents_decomp[i][2]) * 1A  +  ")
    println("$(simpleloopscurrents_decomp[i][3]) * 0.1A")
end

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 1:8
    plot_loop(simpleloops[i], vertex_positions, color = "C0")
end
title("Current: $(simpleloopscurrents[1])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 9:16
    plot_loop(simpleloops[i], vertex_positions, color = "C1")
end
title("Current: $(simpleloopscurrents[9])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 17:20
    plot_loop(simpleloops[i], vertex_positions, color = "C2")
end
title("Current: $(simpleloopscurrents[17])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 21:28
    plot_loop(simpleloops[i], vertex_positions, color = "C3")
end
title("Current: $(simpleloopscurrents[21])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 29:30
    plot_loop(simpleloops[i], vertex_positions, color = "C4")
end
title("Current: $(simpleloopscurrents[29])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 31:38
    plot_loop(simpleloops[i], vertex_positions, color = "C5")
end
title("Current: $(simpleloopscurrents[31])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 39:42
    plot_loop(simpleloops[i], vertex_positions, color = "C6")
end
title("Current: $(simpleloopscurrents[39])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 43:46
    plot_loop(simpleloops[i], vertex_positions, color = "C7")
end
title("Current: $(simpleloopscurrents[43])")

figure(figsize = (6, 6))
plot_vertices(vertex_positions, labels = false)
for i in 47:48
    plot_loop(simpleloops[i], vertex_positions, color = "C8")
end
title("Current: $(simpleloopscurrents[47])")
