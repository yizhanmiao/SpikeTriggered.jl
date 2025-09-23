using Test
using SpikeTriggered

@testset "basics" begin
    @test true
end

@testset "receptive fields" begin
    include("footprint.jl")
end
