#### Binned PSTH
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
    spike_histogram(spike_train::AbstractSpikeTrain{T}, edges::AbstractVector; dtype=Float64, kwargs...) -> Vector

Peri-stimulus time histogram.
`edges` should be length of `n+1` to create a histogram of length `n`.

If `spike_train` is empty, returns `zeros(dtype, n)`.
If `T` is a floating-point type, `dtype` is ignored and `T` is used instead.

## Arguments
- `spike_train`: vector of spike times.
- `edges`: bin edges for the histogram.
- `dtype`: element type of the output (default `Float64`); overridden by `T` when `T <: AbstractFloat`.
- `kwargs...`: forwarded to `FHist.Hist1D`.
"""
function spike_histogram(
    spk::AbstractSpikeTrain{T},
    edges::AbstractVector;
    dtype::Union{Type,Nothing}=Float64,
    kwargs...,
) where {T}
    dtype = T <: AbstractFloat ? T : dtype
    if isempty(spk)
        return zeros(dtype, length(edges) - 1)
    else
        return histogram_fhist(spk, edges; counttype=dtype, kwargs...)
    end
end

# markers can be 2d matrix with shape of [nEdges x nRepeats]
# or specify the repeat dimension
@doc raw"""
    spike_histogram(spike_train::AbstractSpikeTrain, edges::AbstractMatrix; dims=2, kwargs...) -> Matrix

Peri-stimulus time histogram with multiple trials.
`edges` should be shape of `(n+1) × nRepeats` to create a histogram of shape `n × nRepeats`.

## Arguments
- `spike_train`: vector of spike times.
- `edges`: bin-edge matrix; each slice along `dims` is one trial's edges.
- `dims`: dimension along which to slice `edges` (default `2`).
- `kwargs...`: forwarded to the single-trial `spike_histogram`.
"""
function spike_histogram(spk::AbstractSpikeTrain, edges::AbstractMatrix; dims=2, kwargs...)
    return reduce(
        hcat, map(x -> spike_histogram(spk, x; kwargs...), eachslice(edges; dims))
    )
end

@doc raw"""
    spike_histogram(raster::SpikeRaster{T}, edges::AbstractVector; norm=true, kwargs...) -> Vector

PSTH from averaging rasters.

## Arguments
- `raster`: collection of spike trains across trials.
- `edges`: bin edges for the histogram.
- `norm`: if `true` (default), divide counts by the number of trials.
- `kwargs...`: forwarded to the single-train `spike_histogram`.
"""
function spike_histogram(
    raster::SpikeRaster{T}, edges::AbstractVector; norm::Bool=true, kwargs...
) where {T}
    _flatten = reduce(vcat, raster; init=T[])
    _psth = spike_histogram(_flatten, edges; kwargs...)
    return norm ? _psth ./ length(raster) : _psth
end

@doc raw"""
    spike_histogram_center(edges::AbstractVector) -> AbstractVector

Compute the bin centres from a sorted edge vector.

## Arguments
- `edges`: sorted bin edges of length `n+1`.

Returns a vector of length `n`.
"""
spike_histogram_center(edges::AbstractVector) = edges[1:(end - 1)] .+ diff(edges) ./ 2

#### Smoothed PSTH
@doc raw"""
    spike_histogram_smoothed(spike_train::AbstractSpikeTrain, proj::AbstractVector{T}, kernel=gaussian_kernel; norm=true, kwargs...) -> Vector{T}

Generate a smoothed PSTH by convolving spike times with a kernel function.

```math
\text{PSTH}(t) = (h * s)(t) = \sum_i \delta(t - t'_i) h (t - t'_i),\; t \in \text{proj}
```

## Arguments
- `spike_train`: vector of spike times.
- `proj`: time-projection axis (evaluation points).
- `kernel`: smoothing kernel (default `gaussian_kernel`); called as `kernel(Δt; kwargs...)`.
- `norm`: if `true` (default), normalise by the peak absolute value.
- `kwargs...`: forwarded to `kernel`.

!!! note
    This function uses `@floop` internally and benefits from multi-threading.
    Start Julia with multiple threads (e.g. `julia -tauto`) for faster execution.
"""
function spike_histogram_smoothed(
    spk::AbstractSpikeTrain,
    proj::AbstractVector{T},
    kernel::Function=gaussian_kernel;
    norm=true,
    kwargs...,
) where {T}
    isempty(spk) && (return zeros(T, size(proj)))
    _psth = similar(proj, T)
    @floop for idx in eachindex(proj)
        _psth[idx] = sum(kernel.(spk .- proj[idx]; kwargs...))
    end
    return norm ? _psth ./ maximum(abs, _psth) : _psth
end

@doc raw"""
    spike_histogram_smoothed(raster::SpikeRaster{T}, proj::AbstractVector; norm=false, kwargs...) -> Vector

Generate a smoothed PSTH from rasters.

## Arguments
- `raster`: collection of spike trains across trials.
- `proj`: time-projection axis (evaluation points).
- `norm`: if `true`, normalise by peak absolute value; if `false` (default), divide by the number of trials.
- `kwargs...`: forwarded to the single-train `spike_histogram_smoothed`.
"""
function spike_histogram_smoothed(
    raster::SpikeRaster{T}, proj::AbstractVector; norm=false, kwargs...
) where {T}
    _flatten = reduce(vcat, raster; init=T[])
    _psth = spike_histogram_smoothed(_flatten, proj; norm, kwargs...)
    return norm ? _psth : _psth ./ length(raster)
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
