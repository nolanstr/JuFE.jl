#using JuFE
include("/Users/nolanstrauss/bin/JuFE/src/JuFE.jl")
using Test

@testset "JuFE.jl" begin
    
    information = JuFE.ReadYAML("./test.yaml")

    @test information.NELX == 3
    @test information.NX == 3
    @test information.NODAL_COORDS[1].ID == 1
    @test information.NODAL_COORDS[1].x == 0.
    @test information.NODAL_COORDS[1].y == 0.
    @test information.ELEMENTS[1].ID == 1
    @test information.ELEMENTS[1].NODE_IDS == [1, 2]
    @test information.ELEMENTS[1].AREA == 1.
    @test information.ELEMENTS[1].MATERIAL_ID == 1
    @test information.MATERIALS[1].ID == 1
    @test information.MATERIALS[1].E == 10.
    @test information.MATERIALS[1].nu == 0.3

end
