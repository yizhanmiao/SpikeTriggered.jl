using Documenter

if isnothing(get(ENV, "CI", nothing))
    push!(LOAD_PATH, "../src/")

    using Pkg: Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.resolve()
    Pkg.instantiate()
end

using SpikeTriggered

makedocs(;
    modules=[SpikeTriggered],
    sitename="SpikeTriggered.jl",
    #  checkdocs=:exports,
    remotes=nothing,
    format=Documenter.HTML(;
        #   repolink="/dir?ci=SpikeTriggered.jl&name=SpikeTriggered.jl",
        #   disable_git=true,
        #   edit_link=nothing,
        prettyurls=get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Spike Statistics" => "spike_stats.md",
        "Raster and PSTH" => "raster_psth.md",
        "Reverse Correlation" => "reverse_correlation.md",
        #    "Forward Correlation" => "forward_correlation.md",
        "Spatiotemporal Receptive Field" => "strf.md",
        "APIs" => "api.md",
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo="github.com/ZaneMuir/SpikeTriggered.jl.git")
end
