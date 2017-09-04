# Dear Michael,
# thanks a lot again for showing me your very nice setup!
# The cage for GBAR should be 2x2x2 m3 and the region where the gradient should be smaller than 0.02 G/m would be in the middle and would be 50x50x50cm3….the cage should have no walls on the both ends  (the length of the cage could be longer if required). Thank you so much:-)
# Best,
# Paolo

# goal: 2 μT/m

using PyPlot

push!(LOAD_PATH, pwd())
using Coils
using CoilsPlot

g, vertex_positions = cuboid_system([2, 2/5*7, 2], [5, 7, 5],
    skipfaces = [false, false, true, true, false, false])
poi = cuboid_poi([0.6, 0.6, 0.6], [0.0, 0.0, 0.0], [9, 9, 9], filled = false)

# prepare the large cells
front_face_unordered = [ i for i in eachindex(vertex_positions)
    if vertex_positions[i][2] < -1.35 ]
front_face = order_cell(g, front_face_unordered)

back_face_unordered = [ i for i in eachindex(vertex_positions)
    if vertex_positions[i][2] > 1.35 ]
back_face = order_cell(g, back_face_unordered)

initialcells = [front_face, back_face]

elemcurrents = [10.0, 1.0, 0.1]

function design(Bgoal, name)
    simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal,
        verbose = true,
        tolerance = 1e-20,
        bigfloat = false,
        λ = 1e-7, # regularization parameter
        mincurrent = minimum(elemcurrents) / 2,
        initialcells = initialcells)
    simpleloopscurrents_decomp = decompose_currents(simpleloopscurrents, elemcurrents)

    save_report("GBAR_$name", vertex_positions, g, poi, Bgoal, simpleloops,
        simpleloopscurrents, elemcurrents, simpleloopscurrents_decomp,
        levels = linspace(-30, 30, 21))
end

design(x -> [100e-6, 0, 0], "X")
design(x -> [0, 100e-6, 0], "Y")
design(x -> [0, 0, 100e-6], "Z")

design(x -> [100e-6 * x[1], 0, 0], "dBx_dx")
design(x -> [100e-6 * x[2], 0, 0], "dBx_dy")
design(x -> [100e-6 * x[3], 0, 0], "dBx_dz")

design(x -> [0, 100e-6 * x[1], 0], "dBy_dx")
design(x -> [0, 100e-6 * x[2], 0], "dBy_dy")
design(x -> [0, 100e-6 * x[3], 0], "dBy_dz")

design(x -> [0, 0, 100e-6 * x[1]], "dBz_dx")
design(x -> [0, 0, 100e-6 * x[2]], "dBz_dy")
design(x -> [0, 0, 100e-6 * x[3]], "dBz_dz")
