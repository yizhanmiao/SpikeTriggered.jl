# Spike Train Statistics
module SpikeTriggered

using Distributions: cdf, Normal
using DSP: conv
using FFTW: FFTW
using FHist: bincounts, Hist1D
using FLoops: @floop
using GeometryBasics: Point2
using Interpolations: BSpline, interpolate, Linear, Quadratic
using Optim: minimizer, NelderMead, optimize, Options as OptimOptions
using Random: randperm
using SpecialFunctions: erfc, erfcinv
using Statistics: mean, std, var
using StatsBase: weights
using ToeplitzMatrices: Circulant

## general statistics
export AbstractSpikeTrain
export SpikeRaster, AbstractPSTH
export AbstractStimulus, AbstractMarker

export troi2tproj
export modulation_index, dispersion_index
export osi_global, dsi_global

export spike_raster
export spike_histogram, spike_histogram_center
export spike_histogram_smoothed
export spike_spectrum_power
public gaussian_kernel

## reverse correlation
export spike_triggered_average
export spike_triggered_average_rectified
export spike_triggered_average_zscore
public spike_triggered_average_bootstrap

## forward correlation
export spike_raster_groupby

## bursts
export spike_detect_burst
export spike_detect_burst_trn, spike_detect_burst_lgn
export spike_split_burst, spike_split_burst_cardinal
export spike_burst_pattern

## strfs
export AbstractSTRF, AbstractTRF, AbstractSRF
export strf_, strf_hstack, strf_vstack

# srf
export GaussianEllipse, srf_gaussian_fit
export srf_schiller_overlap_index
export srf_footprint, srf_footprint_mask
export srf_pixel_cloud

public area, make_heatmap, make_ellipse
public floodfill!, image2dataset
public collapse_footprint
public benjamini_hochberg_constant

# trf
export trf_segment_biphasic
export trf_segment_statistics
#TODO: export trf_isbiphasic
export trf_polarity_score

const SpikeRaster{T} = AbstractVector{<:AbstractVector{T}}
const AbstractSpikeTrain{T} = AbstractVector{T}
const AbstractPSTH{T} = AbstractVecOrMat{T}    # [nTimepoints] or [nTimepoints x nRepeats]
const AbstractStimulus{T} = AbstractMatrix{T}  # [nDimensions x nTimepoints]
const AbstractMarker{T} = AbstractVecOrMat{T}  # [nTimepoints] or [nTimepoints x nRepeats]

include("statistics.jl")
include("raster.jl")
include("psth.jl")
include("spike_triggered_average.jl") # a.k.a. spike-triggered averaging
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

# Deprecated names
@deprecate burst_detect spike_detect_burst
@deprecate burst_detect_trn spike_detect_burst_trn
@deprecate burst_detect_lgn spike_detect_burst_lgn
@deprecate split_tonic_burst spike_split_burst
@deprecate split_tonic_cardinal spike_split_burst_cardinal
@deprecate burst_interpolate spike_burst_pattern
@deprecate orientation_selectivity_index_global osi_global
@deprecate direction_selectivity_index_global dsi_global
@deprecate event_triggered_raster spike_raster_groupby
@deprecate get_histogram_center spike_histogram_center
@deprecate spike_triggered_average_suite spike_triggered_average_rectified
@deprecate get_trf_polarity_score trf_polarity_score
@deprecate get_footprint_map srf_footprint
@deprecate get_footprint_mask srf_footprint_mask
@deprecate get_pixel_cloud srf_pixel_cloud

end
