using Test
using SpikeTriggered
using Aqua

@testset "Aqua" begin
    Aqua.test_all(SpikeTriggered)
end

@testset "basics" begin
    @test true
end

@testset "receptive fields" begin
    include("footprint.jl")
end

@testset "spike statistics" begin
    include("psth.jl")
    include("burst.jl")
end