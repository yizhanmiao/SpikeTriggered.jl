export spike_histogram, get_histogram_center
export spike_histogram_smoothed

#### Binned PSTH

@doc raw"""
    histogram_gsl(u_arr, edges) where {T <: Real}

Count histogram by providing the edges. If provides `n+1` edges,
it would return `n` length histogram.
Edges would work as: left bound <= value < right bound.

# Arguments
- `u_arr`: a vector of values
- `edges`: a vector of edges

Note: StatsBase.Histogram is slower than GSL but with less memory footprint.
"""
function histogram_gsl(u_arr::AbstractVector{T}, edges) where {T<:Real}
    u_arr = Cdouble.(u_arr)
    edges = Cdouble.(edges)

    edges = sort(edges)

    n = length(edges) - 1
    gsl_hist = GSL.histogram_alloc(n)
    GSL.histogram_set_ranges(gsl_hist, edges, n + 1)

    for idx in eachindex(u_arr)
        @inbounds GSL.histogram_increment(gsl_hist, u_arr[idx])
    end

    myhist = zeros(Int, n)
    for idx in 1:n
        @inbounds myhist[idx] = GSL.histogram_get(gsl_hist, idx - 1)
    end

    GSL.histogram_free(gsl_hist)

    return myhist
end

@doc raw"""
    histogram_fhist(u_arr, edges; kwargs...) -> Vector

Count histogram using `FHist.Hist1D`, which has similar performance as `GSL`.

`kwargs` will be relayed to `FHist.Hist1D`;
for example, one could set element type of the histogram to `Float32`
by passing `counttype=Float32`.
"""
function histogram_fhist(u_arr::AbstractVector, edges; kwargs...)
    h = Hist1D(u_arr; binedges=edges, kwargs...)
    return bincounts(h)
end

@doc raw"""
    spike_histogram(spike_train::AbstractVector{T}, edges::AbstractVector; dtype::Type, kwargs...) -> Vector{T}

Peri-stimulus time histogram.
`edges` should be length of `n+1` to create a histogram of length `n`.

If `spike_train` is empty, it will return `zeros(T, n)`

If `counttype` is not specified, the returned vector will have same type `T`
"""
function spike_histogram(
    spk::AbstractSpikeTrain{T},
    edges::AbstractVector;
    dtype::Union{Type,Nothing}=nothing,
    kwargs...,
) where {T<:Real}
    dtype = isnothing(dtype) ? T : dtype
    if isempty(spk)
        zeros(dtype, length(edges) - 1)
    else
        histogram_fhist(spk, edges; counttype=dtype, kwargs...)
    end
end

# markers can be 2d matrix with shape of [nEdges x nRepeats]
# or specify the repeat dimension
@doc raw"""
    spike_histogram(spike_train, edges::AbstractMatrix; kwargs...) -> Matrix{T}

Peri-stimulus time histogram with multiple trials.
`edges` should be shape of [`n+1` x nRepeats] to create a histogram of shape [`n` x nRepeats].

If `dims` is specified, trial dimension will be overwrite.

Additional `kwargs` will be passed to `spike_histogram`.
"""
function spike_histogram(spk::AbstractSpikeTrain, edges::AbstractMatrix; dims=2, kwargs...)
    return reduce(
        hcat, map(x -> spike_histogram(spk, x; kwargs...), eachslice(edges; dims))
    )
end

@doc raw"""
    spike_histogram(raster::SpikeRaster{T}, args...; norm=true, kwargs...) -> VecOrMat{T}

PSTH from averaging rasters.

If `norm` is `true`, PSTH will be normalized by the length of trials.
"""
function spike_histogram(
    raster::SpikeRaster{T}, args...; norm::Bool=true, kwargs...
) where {T<:Real}
    _flatten = reduce(vcat, raster; init=T[])
    _psth = spike_histogram(_flatten, args...; kwargs...)
    return norm ? _psth ./ length(raster) : _psth
end

@doc raw"""
    get_histogram_center(edges)

get the center of edge vector. make sure edges are sorted.
"""
get_histogram_center(edges::AbstractVector) = edges[1:(end - 1)] .+ diff(edges) ./ 2

#### Smoothed PSTH

@deprecate spike_filter(spk, proj, kernel; norm_by=nothing, kwargs...) spike_histogram_smoothed(
    spk, proj, kernel; norm=isnothing(norm_by), kwargs...
) false
@deprecate spike_filter_gaussian(spk_or_raster, proj; kwargs...) spike_histogram_smoothed(
    spk_or_raster, proj; kwargs...
) false

@doc raw"""
    spike_histogram_smoothed(spike_train::SpikeTrain{T}, projection::AbstractArray, kernel::Function=gaussian_kernel; norm=true, kwargs...) -> Vector{T}

Generating smoothed curve from spike trains. Equivalent to convolution.

```math
\text{PSTH}(t) = (h * s)(t) = \sum_i \delta(t - t'_i) h (t - t'_i),\; t \in \text{proj}
```

If `norm` is `true`, the results will be normalized by the maximum value.

All `kwargs` will be passed to `kernel` function.
"""
function spike_histogram_smoothed(
    spk::AbstractSpikeTrain{T},
    proj::AbstractArray,
    kernel::Function=gaussian_kernel;
    norm=true,
    kwargs...,
) where {T<:Real}
    isempty(spk) && (return zeros(T, size(proj)))
    _psth = similar(proj, T)
    @floop for idx in eachindex(proj)
        _psth[idx] = sum(kernel.(spk .- proj[idx]; kwargs...))
    end
    return norm ? _psth ./ maximum(abs, _psth) : _psth
end

@doc raw"""
    spike_histogram_smoothed(raster, args...; kwargs...) where {T <: Real} -> Vector{T}

Generating smoothed PSTH from rasters.

If `norm` is `true`, PSTH will be normalized to the maximum number;
otherwise, it will be divided by the number of trials in the raster.
"""
function spike_histogram_smoothed(
    raster::SpikeRaster{T}, args...; norm=false, kwargs...
) where {T<:Real}
    _flatten = reduce(vcat, raster; init=T[])
    _N = length(raster)
    if norm
        spike_histogram_smoothed(_flatten, args...; norm=true, kwargs...)
    else
        spike_histogram_smoothed(_flatten, args...; norm=false, kwargs...) ./ _N
    end
end

@doc raw"""
    gaussian_kernel(x::T; σ::T=0.005)

Inline [Gaussian function](https://en.wikipedia.org/wiki/Gaussian_filter).

```math
g(x) = \frac{\exp(-\frac{x^2}{2 σ^2})}{\sqrt{2\pi} σ}
```
"""
@inline gaussian_kernel(x::T; σ::T=0.005) where {T<:AbstractFloat} =
    exp(x^2 / (-2 * σ^2)) / (σ * sqrt(2 * T(π)))
