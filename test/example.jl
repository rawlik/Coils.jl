# code from the example notebook for automatic testing

using PyPlot

using Coils
using Coils.CoilsPlot

g, vertex_positions = cuboid_system([1, 1, 1], [3, 3, 3])

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions)
plot_edges(ax, g, vertex_positions)

poi = cuboid_poi([0.75, 0.75, 0.75], [0.0, 0.0, 0.0], [10, 10, 10], filled = false)

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_system(ax, g, vertex_positions, poi)

Bgoal(x) = [0, 100e-6, 0]

cells = find_cells(g)

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
plot_cells(ax, cells, vertex_positions)

M = system_matrix(poi, vertex_positions, cells)

Bpoi = vcat(Bgoal.(poi)...)

optI = M \ Bpoi

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
current_norm = 1000 / maximum(abs, optI)
plot_cells(ax, cells, vertex_positions, optI * current_norm)

edgecurrents = find_edgecurrents(g, cells, optI)

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
plot_edge_currents(ax, g, edgecurrents, vertex_positions)

simpleloops, simpleloopscurrents = find_all_simpleloops(g, edgecurrents)

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
plot_loop(ax, simpleloops[1], vertex_positions, color = "orange")
plot_loop(ax, simpleloops[2], vertex_positions, color = "orange")
title("Current: $(simpleloopscurrents[1] * current_norm)")

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
plot_loop(ax, simpleloops[3], vertex_positions, color = "green")
plot_loop(ax, simpleloops[4], vertex_positions, color = "green")
title("Current: $(simpleloopscurrents[3] * current_norm)")

fig, ax = subplots(figsize = (6, 6))
plot_deviation_histogram(ax, poi, simpleloops, simpleloopscurrents, vertex_positions, Bgoal)

plot_loops_field(cells, optI, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-0.5, 0.5], Bref = Bgoal, levels = -10:10)

g, vertex_positions = cuboid_system([1, 2, 1], [3, 6, 3],
    skipfaces = [false, false, true, true, false, false])
poi = cuboid_poi([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [10, 10, 10], filled = true)

# the wires used to construct the coil
elemcurrents = [10.0, 1.0, 0.1]

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_system(ax, g, vertex_positions, poi)

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

fig, ax = subplots(figsize = (6, 6))
plot_deviation_histogram(ax, poi, simpleloops, simpleloopscurrents_real, vertex_positions, Bgoal)

plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-1, 1], Bref = Bgoal, levels = -20:4:20)

plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions, [1, 2, 3], 0;
        n = 50, spanA = [-0.5, 0.5], spanB = [-1, 1], levels = -50:5:50)

for i in eachindex(simpleloopscurrents)
    print("#$(lpad(i,2)): ")
    print("$(lpad(simpleloopscurrents[i], 20)) = ")
    print("$(lpad(simpleloopscurrents_decomp[i][1], 2)) * 10A  +  ")
    print("$(simpleloopscurrents_decomp[i][2]) * 1A  +  ")
    println("$(simpleloopscurrents_decomp[i][3]) * 0.1A")
end

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
for i in 1:8
    plot_loop(ax, simpleloops[i], vertex_positions, color = "C0")
end
title("Current: $(simpleloopscurrents[1])")

fig = figure(figsize = (6, 6))
ax = fig.add_subplot(projection="3d")
plot_vertices(ax, vertex_positions, labels = false)
for i in 9:16
    plot_loop(ax, simpleloops[i], vertex_positions, color = "C1")
end
title("Current: $(simpleloopscurrents[9])")
