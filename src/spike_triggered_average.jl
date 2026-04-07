@doc raw"""
    spike_triggered_average(stimulus, psth; n=10)
    spike_triggered_average(stimulus, spike_train, marker; n=10)
    spike_triggered_average(stimulus, raster, marker; n=10)

Compute spike-triggered average (STA), which can be reshaped into a 2D/3D matrix
with `make_strf` and `hstack_strf`.

If `psth` is a 2D matrix (nTimepoints × nRepeats),
the average of PSTH across trials will be used.

```math
A = \frac{1}{N} \sum_{n=1}^{N}\vec{s}(t_n)
```

## Arguments
- `stimulus::AbstractStimulus`: stimulus matrix (dimensions × timepoints).
- `psth`: PSTH vector or matrix, or a spike train / raster with `marker`.
- `marker`: stimulus onset marker (used with spike train / raster).
- `n`: number of preceding frames to include (default `10`).

Returns a flattened vector of length `D × n`.

!!! tip
    For a high-contrast STA, rectify the stimulus before computing
    (e.g. `clamp.(stimulus, 0, Inf)` for ON, `clamp.(stimulus, -Inf, 0)` for OFF).
    See [`spike_triggered_average_rectified`](@ref) for a convenience wrapper.
"""
function spike_triggered_average(
    stimulus::AbstractStimulus, psth::AbstractVector; n::Integer=10
)
    if size(stimulus, 2) != length(psth)
        throw(
            ArgumentError(
                "stimulus timepoints ($(size(stimulus, 2))) must match PSTH length ($(length(psth)))",
            ),
        )
    end

    resp = Circulant(psth)
    sta = stimulus * resp[:, [1; end:-1:(end - n + 2)]] ./ sum(psth)  # [Dimensions x N]
    return sta[:]
end

function spike_triggered_average(
    stimulus::AbstractStimulus, psth::AbstractMatrix; kwargs...
)
    return spike_triggered_average(stimulus, mean(psth; dims=2)[:]; kwargs...)
end

function spike_triggered_average(
    stimulus::AbstractStimulus,
    spike_or_raster::Union{AbstractSpikeTrain,SpikeRaster},
    marker::AbstractMarker;
    kwargs...,
)
    return spike_triggered_average(
        stimulus, spike_histogram(spike_or_raster, marker); kwargs...
    )
end

## STA diff-on-off

function _rectify_stimulus(stim; bias=:diff)
    if bias == :on
        return clamp.(stim, 0, Inf)
    elseif bias == :off
        return clamp.(stim, -Inf, 0)
    else
        return copy(stim)
    end
end

@doc raw"""
    spike_triggered_average_rectified(stimulus, args...; kwargs...) -> NamedTuple{(:diff, :on, :off)}

Compute STAs for the full stimulus (diff), ON-rectified, and OFF-rectified versions.

The ON map zeroes out negative stimulus values; the OFF map zeroes out positive values.

## Arguments
- `stimulus::AbstractStimulus`: stimulus matrix (dimensions × timepoints).
- `args...`, `kwargs...`: forwarded to `spike_triggered_average`.

Returns a named tuple `(; difference, bright, dark)`.
"""
function spike_triggered_average_rectified(stimulus, args...; kwargs...)
    return (;
        difference=spike_triggered_average(stimulus, args...; kwargs...),
        bright=spike_triggered_average(
            _rectify_stimulus(stimulus; bias=:on), args...; kwargs...
        ),
        dark=spike_triggered_average(
            _rectify_stimulus(stimulus; bias=:off), args...; kwargs...
        ),
    )
end

## STA bootstrapping

@doc raw"""
    spike_triggered_average_bootstrap(stimulus, psth; n=10, bootstrap=0) -> Matrix{Float64}

Bootstrap the STA by circularly shifting the PSTH relative to the stimulus.

When `bootstrap=0` (default), all possible shifts are used (full permutation test).
Otherwise, `bootstrap` random shift offsets are sampled.

## Arguments
- `stimulus::AbstractStimulus`: stimulus matrix (dimensions × timepoints).
- `psth::AbstractPSTH`: PSTH vector or matrix.
- `n`: number of preceding frames (default `10`).
- `bootstrap`: number of bootstrap iterations; `0` for exhaustive (default `0`).

Returns a `(D×n) × nIterations` matrix of bootstrapped STAs.

!!! warning
    Prefer the spike-train method for final analysis — circular shifting
    of a PSTH can introduce artefacts.
"""
function spike_triggered_average_bootstrap(
    stimulus::AbstractStimulus, psth::AbstractPSTH; n::Integer=10, bootstrap::Integer=0
)
    @warn "please use bootstrapping with spike times for final analysis."

    (D, N) = size(stimulus)
    _bootstrap_offset = if bootstrap == 0  ## full range
        collect(range(n + 1; step=1, length=N - n))
    else
        randperm(N - n)[1:min(N, bootstrap)] .+ n
    end

    bsta = Matrix{Float64}(undef, D * n, length(_bootstrap_offset))
    @floop for (idx, offset) in enumerate(_bootstrap_offset)
        @inbounds bsta[:, idx] .= spike_triggered_average(
            stimulus, circshift(psth, -offset); n
        )
    end
    return bsta
end

@doc raw"""
    spike_triggered_average_bootstrap(stimulus, spike_train, marker; n=10, bootstrap=1000, skip_n=false) -> Matrix{Float64}

Bootstrap the STA by randomly shifting spike times within the stimulus duration.

## Arguments
- `stimulus::AbstractStimulus`: stimulus matrix (dimensions × timepoints).
- `spike_train::AbstractSpikeTrain`: vector of spike times.
- `marker::AbstractMarker`: stimulus onset times.
- `n`: number of preceding frames (default `10`).
- `bootstrap`: number of bootstrap iterations (default `1000`).
- `skip_n`: if `true`, exclude shifts shorter than `n` average frames (default `false`).

Returns a `(D×n) × bootstrap` matrix of bootstrapped STAs.
"""
function spike_triggered_average_bootstrap(
    stimulus::AbstractStimulus,
    spike_train::AbstractSpikeTrain,
    marker::AbstractMarker;
    n::Integer=10,
    bootstrap::Integer=1000,
    skip_n::Bool=false,
)
    (D, _N) = size(stimulus)
    if bootstrap <= 0
        throw(ArgumentError("bootstrap must be a positive integer, got $bootstrap"))
    end

    _lo, _hi = extrema(marker)
    _duration = _hi - _lo
    _t_shift = if skip_n
        _delta_t = mean(diff(marker[:]))
        rand(bootstrap) .* (_duration - _delta_t * n) .+ _delta_t * n
    else
        rand(bootstrap) .* _duration
    end

    bsta = Matrix{Float64}(undef, D * n, bootstrap)
    @floop for (idx, offset) in enumerate(_t_shift)
        @inbounds bsta[:, idx] .= spike_triggered_average(
            stimulus, ((spike_train .+ offset) .% _duration .+ marker[1]), marker; n
        )
    end
    return bsta
end

@doc raw"""
    spike_triggered_average_zscore(args...; n=10, bootstrap=1000) -> Vector

Estimate z-scores of the STA against a bootstrap null distribution.

Computes the STA and a bootstrapped distribution, then returns the
element-wise z-score: ``(STA - \mu) / \sigma``.

## Arguments
- `args...`: forwarded to both `spike_triggered_average` and `spike_triggered_average_bootstrap`.
- `n`: number of preceding frames (default `10`).
- `bootstrap`: number of bootstrap iterations (default `1000`).

Returns a vector of z-scores with length `D × n`.
"""
function spike_triggered_average_zscore(args...; n::Integer=10, bootstrap::Integer=1000)
    sta = spike_triggered_average(args...; n)
    bsta = spike_triggered_average_bootstrap(args...; n, bootstrap)

    μ = mean(bsta; dims=2)[:]
    σ = std(bsta; dims=2)[:]

    return (sta .- μ) ./ σ
end
