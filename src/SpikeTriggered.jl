# Spike Train Statistics
module SpikeTriggered

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
using DSP: conv
using Optim: optimize, NelderMead, Options as OptimOptions, minimizer
using Interpolations: interpolate, BSpline, Linear, Quadratic

const SpikeRaster{T} = AbstractVector{<:AbstractVector{T}}
const AbstractSpikeTrain{T} = AbstractVector{T}
const AbstractPSTH{T} = AbstractVecOrMat{T}    # [nTimepoints] or [nTimepoints x nRepeats]
const AbstractStimulus{T} = AbstractMatrix{T}  # [nDimensions x nTimepoints]
const AbstractMarker{T} = AbstractVecOrMat{T}  # [nTimepoints] or [nTimepoints x nRepeats]

export SpikeRaster, AbstractSpikeTrain

include("statistics.jl")
include("raster.jl")
include("psth.jl")
include("spike_triggered_average.jl") # reverse correlation
include("forward_correlation.jl")
include("burst.jl")
include("spectrum.jl")
include("drifting_grating.jl")
include("strf/strf.jl")

#TODO: review the functions below and make documentations
include("jpsth.jl")
include("reliability.jl")
include("correlogram.jl")
include("spike_triggered_covariance.jl")

# additional waveform related functions
include("Waveforms/Waveforms.jl")

# using Reexport: @reexport
# @reexport using .Waveforms

end
