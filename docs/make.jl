using Documenter
using Literate
using SpikeTriggered

EXAMPLE = joinpath(@__DIR__, "src", "examples")

Literate.markdown.([joinpath(EXAMPLE, "raster_psth.jl")], EXAMPLE; documenter=true)

# include("changelog.jl")

include("pages.jl")

makedocs(;
    modules=[SpikeTriggered],
    sitename="SpikeTriggered.jl",
    remotes=nothing,
    checkdocs=:exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        collapselevel=1,
        repolink="https://taro-station.usc.edu:20443/AnalysisPipeline/SpikeTriggered.jl",
    ),
    pages=pages,
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo="github.com/yizhanmiao/SpikeTriggered.jl.git")
else
    @info "You are build docs locally."
end
