using PyPlot

push!(LOAD_PATH, pwd())
using Coils
using CoilsPlot

g, vertex_positions = cuboid_system([0.4, 0.5, 0.4], [4, 5, 4],
    skipfaces = [false, false, true, false, false, false])
poi = cuboid_poi([0.3, 0.3, 0.3], [0.0, (0.5 - 0.3 - 0.1) / 2, 0.0], [9, 9, 9], filled = true)

function design(Bgoal, name)
    simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal,
        verbose = true, tolerance = 1e-20, bigfloat = true,
        initialcells = [[1, 2, 3, 4, 5, 32, 47, 62, 80, 79, 78, 77, 76, 61, 46, 31, 1]])

    elemcurrents = [1.0, 0.1, 0.01]
    simpleloopscurrents_decomp = decompose_currents(simpleloopscurrents, elemcurrents)

    save_report("small_$name", vertex_positions, g, poi, Bgoal, simpleloops,
        simpleloopscurrents, elemcurrents, simpleloopscurrents_decomp)
end

design(x -> [100e-6, 0, 0], "X")
design(x -> [0, 100e-6, 0], "Y")
design(x -> [0, 0, 100e-6], "Z")

design(x -> [10e-6 * x[1], 0, 0], "dBx_dx")
design(x -> [10e-6 * x[2], 0, 0], "dBx_dy")
design(x -> [10e-6 * x[3], 0, 0], "dBx_dz")
