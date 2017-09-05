import Base.Test.@testset
import Base.Test.@test

using LightGraphs
using PyPlot
using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport pySFC

push!(LOAD_PATH, pwd())
using Coils


const tol = 1e-10


@testset "Biot-Savart tests" begin
    loop_size = [1, 1]
    current = 100
    centre = [0, 0, 0]
    p = [[0,0,0], [0, 0, 20], [0.2, 0.4, 0.5], [-0.4, -0.4, -0.3],
        [ -0.691412, -0.568035, -0.404835], [-0.38749, -0.737193, 0.295644],
        [-2, -2, -2], [2, 2, 2], [2, 2, 0]]

    # python coil
    pycoil = pySFC.Coil(center = centre, current = current, size_x = loop_size[1], size_y = loop_size[2],
        rotation = eye(3))
    pyfield = [ pycoil[:B](p...) for p in p ]

    m = [loop_size / 2; 0]
    vertex_positions = [[1,1,0] .* m, [1,-1,0] .* m, [-1,-1,0] .* m, [-1,1,0] .* m]
    # g = DiGraph()
    # add_vertices!(g, length(vertex_positions))
    # add_edge!(g, 1, 2)
    # add_edge!(g, 2, 3)
    # add_edge!(g, 3, 4)
    # add_edge!(g, 4, 1)
    juliafield = [ biotsavart_cell(p, vertex_positions, [1, 2, 3, 4, 1]) * current for p in p ]
    # juliafield = [ field_loops(p, [[1, 2, 3, 4, 1]], [current], vertex_positions) for p in p ]
    # Bgoal(x) = [0, 0, 50]
    # simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal)

    # valid only for p on axis of the coil in a large distance from the coil
    analfield = [ μ0 * current * prod(loop_size) / (2π * p[3]^3) for p in p ]

    # test in the field in the loop centre
    @test all(abs.(juliafield[1] .- pyfield[1]) .< tol)

    # test in the field on the axis
    @test all(abs.(juliafield[2] .- analfield[2]) .< 1e-8)
    @test all(abs.(juliafield[2] .- pyfield[2]) .< tol)

    # test in the field in an arbitrary point
    @test all(abs.(juliafield[3] .- pyfield[3]) .< tol)
    @test all(abs.(juliafield[4] .- pyfield[4]) .< tol)
    @test all(abs.(juliafield[5] .- pyfield[5]) .< tol)
    @test all(abs.(juliafield[6] .- pyfield[6]) .< tol)
    @test all(abs.(juliafield[7] .- pyfield[7]) .< tol)

    # @test begin
    # @test foo(zeros(2)) == 4
    # @test foo(ones(4)) == 15

    g = linspace(-2,2,8)
    function testpoint(p)
        pyfield = pycoil[:B_response](p...) * current
        juliafield = biotsavart_cell(p, vertex_positions, [1, 2, 3, 4, 1]) * current
        abs.(juliafield .- pyfield) .< 1e-15
    end

    for p in [[x,y,z] for x in g for y in g for z in g]
        @test testpoint(p)
    end
end


@testset "Solving a system" begin
    vertex_positions = [[1,1,0], [1,-1,0], [-1,-1,0], [-1,1,0]]
    g = DiGraph()
    add_vertices!(g, length(vertex_positions))
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 1)
    poi = [[0,0,0]]
    Bgoal(x) = [0, 0, 50]
    simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal)

    field = field_loops(poi[1], simpleloops, simpleloopscurrents, vertex_positions)
    @test abs(field[3] - 50) < tol
end
