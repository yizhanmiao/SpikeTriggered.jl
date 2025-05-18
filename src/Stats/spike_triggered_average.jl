export spike_triggered_average
export spike_triggered_average_suite
export spike_triggered_average_zscore

@doc raw"""
    spike_triggered_average(stimulus, psth; n=10)
    spike_triggered_average(stimulus, spike_train, marker; n=10)
    spike_triggered_average(stimulus, raster, marker; n=10)

Get the STA, which can be reshaped into 2d/3d matrix
with `make_strf` and `hstack_strf`.

If `psth` is 2d matrix (nTimepoints x nRepeats),
the average of PSTH across trials will be used.

```math
A = \frac{1}{N} \sum_{n=1}^{N}\vec{s}(t_n)
```
"""
function spike_triggered_average(
    stimulus::AbstractStimulus, psth::AbstractVector; n::Integer=10
)
    @assert size(stimulus, 2) == length(psth) "stimulus length must matches with PSTH size."

    resp = Circulant(psth)
    sta = stimulus * resp[:, [1; end:-1:(end - n + 2)]] ./ sum(psth)  # [Dimensions x N]
    return sta[:]
end

function spike_triggered_average(
    stimulus::AbstractStimulus, psth::AbstractMatrix; kwargs...
)
    spike_triggered_average(stimulus, mean(psth; dims=2)[:]; kwargs...)
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

function _simple_stimulus_nonlinearity(stim; bias=:diff)
    video = copy(stim)
    if bias == :on
        video[video .< 0] .= 0
    elseif bias == :off
        video[video .> 0] .= 0
    end
    video
end

@doc raw"""
    spike_triggered_average_suite(args...; kwargs...) -> (; diff, on, off)

Create STAs of difference map, on map and off map; returns a named tuple.
"""
function spike_triggered_average_suite(stimulus, args...; kwargs...)
    (;
        diff=spike_triggered_average(stimulus, args...; kwargs...),
        on=spike_triggered_average(
            _simple_stimulus_nonlinearity(stimulus; bias=:on), args...; kwargs...
        ),
        off=spike_triggered_average(
            _simple_stimulus_nonlinearity(stimulus; bias=:off), args...; kwargs...
        ),
    )
end

## STA bootstrapping

function spike_triggered_average_bootstrap(
    stimulus::AbstractStimulus{T}, psth::AbstractPSTH; n::Integer=10, bootstrap=0
) where {T}
    @warn "please use bootstrapping with spike times for final analysis."

    (D, N) = size(stimulus)
    _bootstrap_offset = if bootstrap == 0  ## full range
        range(n+1; step=1, length=N-n)
    else
        randperm(N-n)[1:min(N, bootstrap)] .+ n
    end

    bsta = Matrix{T}(undef, D*n, length(_bootstrap_offset))
    @floop for (idx, offset) in enumerate(_bootstrap_offset)
        @inbounds bsta[:, idx] .= spike_triggered_average(
            stimulus, circshift(psth, -offset); n
        )
    end
    bsta
end

function spike_triggered_average_bootstrap(
    stimulus::AbstractMatrix{T},
    spike_train::AbstractSpikeTrain,
    marker::AbstractMarker;
    n::Integer=10,
    bootstrap=1000,
    skip_n=false,
) where {T}
    (D, _N) = size(stimulus)
    @assert bootstrap > 0 "bootstrap must be positive value"

    _duration = abs(-(extrema(marker)...))
    _t_shift = if skip_n
        _delta_t = mean(diff(marker[:]))
        rand(bootstrap) .* (_duration - _delta_t * n) .+ _delta_t * n
    else
        rand(bootstrap) .* _duration
    end

    bsta = Matrix{T}(undef, D*n, bootstrap)
    @floop for (idx, offset) in enumerate(_t_shift)
        @inbounds bsta[:, idx] .= spike_triggered_average(
            stimulus, ((spike_train .+ offset) .% _duration .+ marker[1]), marker; n
        )
    end
    bsta
end

@doc raw"""
    spike_triggered_average_zscore(args...; kwargs...)

Estimate z-scores of STA.

Same `args...` as `spike_triggered_average`.

If `n` is specified, only `n` frames are included for zscore (same as sta, default as `10``).

`bootstrap` set the bootstrapping iterations, default as `1000`.
"""
function spike_triggered_average_zscore(args...; n=10, bootstrap=1000)
    sta = spike_triggered_average(args...; n)
    bsta = spike_triggered_average_bootstrap(args...; n, bootstrap)

    μ = mean(bsta; dims=2)[:]
    σ = std(bsta; dims=2)[:]

    (sta .- μ) ./ σ
end

### CUDA version
function cu_spike_triggered_average_bootstrap(
    stimulus::AbstractMatrix, psth; n::Integer=10, bootstrap=0
)
    #TODO
end

function cu_spike_triggered_average_zscore(stimulus, psth; n=10, bootstrap=1)
    #TODO
end
