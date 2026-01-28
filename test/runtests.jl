using Test
using SpikeTriggered

@testset "basics" begin
    @test true
end

@testset "receptive fields" begin
    include("footprint.jl")
end

@testset "spike statistics" begin
    include("psth.jl")
end