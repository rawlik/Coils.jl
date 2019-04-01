import Test.@testset
import Test.@test

using LightGraphs

push!(LOAD_PATH, pwd())
using Coils


const tol = 1e-10

@testset "Solving a system" begin
    g, vertex_positions = cuboid_system([1, 1, 1], [2, 2, 2])
    poi = [[0,0,0]]
    Bgoal(x) = [0, 0, 50]
    simpleloops, simpleloopscurrents = solve_system(g, vertex_positions, poi, Bgoal)

    field = field_loops(poi[1], simpleloops, simpleloopscurrents, vertex_positions)
    @test abs(field[3] - 50) < tol
end

@testset "Example code with plotting" begin
    include("example.jl")
end
