# Spike Train Statistics
module Stats

using GSL: GSL
using FFTW: FFTW
using Random: randperm
using Statistics: mean, var, std
using StatsBase: weights
using SpecialFunctions: erfc, erfcinv
using Distributions: cdf, Normal
using FLoops: @floop
using FHist: Hist1D, bincounts
using ToeplitzMatrices: Circulant
using GeometryBasics: Point2
using Optim: optimize, NelderMead, Options as OptimOptions

const SpikeRaster{T} = Vector{Vector{T}}
const AbstractSpikeTrain{T} = AbstractVector{T}
const AbstractPSTH{T} = AbstractVecOrMat{T}    # [nTimepoints] or [nTimepoints x nRepeats]
const AbstractStimulus{T} = AbstractMatrix{T}  # [nDimensions x nTimepoints]
const AbstractMarker{T} = AbstractVecOrMat{T}  # [nTimepoints] or [nTimepoints x nRepeats]

export SpikeRaster, AbstractSpikeTrain
export spike_triggered_average,
    spike_triggered_average_zscore, spike_triggered_average_suite
export srf_gaussian_fit, Gaussian2DSimplex, Gaussian2DSimplexInit

include("misc.jl")
include("raster.jl")
include("psth.jl")
include("spike_triggered_average.jl")
include("footprint.jl")
include("burst.jl")
include("spectrum.jl")
include("strf/strf.jl")

#TODO: review the functions below and make documentations
include("jpsth.jl")
include("reliability.jl")
include("correlogram.jl")
include("spike_triggered_covariance.jl")

end
