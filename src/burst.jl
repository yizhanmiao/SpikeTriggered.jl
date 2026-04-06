@doc raw"""
    spike_detect_burst(spk::AbstractSpikeTrain; t_silence=(0.07, 0), t_isi=0.03, nofs=nothing, keep_index=false) -> Vector{Vector{T}}

Burst firing identification based on the rules:

- a burst would followed/preceded by a silence `t_silence` time;
if `t_silence` is a 2 element tuple, the `t_silence[1]` will set the preceded time and
`t_silence[2]` will be the followed time.
- each spike spaced at most `t_isi` time, otherwise terminated
- (optional) at least `nofs[1]` spikes within the first `nofs[2]` time

If `keep_index` is `true`, will return indices of spikes within the input;
otherwise, the spike times will be returned.

References:
- [Vaingankar et al 2012](https://doi.org/10.3389/fnint.2012.00118)
"""
function spike_detect_burst(
    spk::AbstractSpikeTrain{T};
    t_silence::Union{Real,NTuple{2,Real}}=(0.07, 0),
    t_isi::Real=0.03,
    nofs::Union{Nothing,Tuple{Integer,Real}}=nothing,
    keep_index=false,
) where {T<:Real}

    # initialization
    t_pre, t_post = if t_silence isa Tuple
        t_silence
    else
        (t_silence, t_silence)
    end
    burst_list = Vector{Vector{T}}()
    burst_index_list = Vector{Vector{Int64}}()
    current_burst = Vector{T}()
    current_burst_index = Vector{Int64}()
    flag_in_burst = false

    # given the spike numbers are usually large, it is okay to skip the first and the last bursts, if any.
    for idx in 2:(length(spk) - 1)
        if flag_in_burst
            # continue
            if spk[idx + 1] - spk[idx] <= t_isi
                push!(current_burst, spk[idx])
                push!(current_burst_index, idx)

                # the last one
            else
                push!(current_burst, spk[idx])
                push!(current_burst_index, idx)
                flag_in_burst = false
                # check the post silence condition
                if spk[idx + 1] - spk[idx] >= t_post
                    push!(burst_list, current_burst)
                    push!(burst_index_list, current_burst_index)
                end
            end

            # seek for potential burst start
        else
            if spk[idx] - spk[idx - 1] >= t_pre && spk[idx + 1] - spk[idx] <= t_isi
                flag_in_burst = true
                current_burst = Vector{T}([spk[idx]])
                current_burst_index = Vector{Int64}([idx])
            end
        end
    end

    rez = if !isnothing(nofs)
        _mask = map(burst_list) do burst
            length(burst) >= nofs[1] && sum(diff(burst[1:nofs[1]])) <= nofs[2]
        end
        burst_list[_mask], burst_index_list[_mask]
    else
        burst_list, burst_index_list
    end

    return keep_index ? rez[2] : rez[1]
end

@doc raw"""
    spike_detect_burst_trn(spk; kwargs...)

≥5 spikes within 70 ms, spaced ≤30 ms apart following ≥70 ms
of silence; bursts were terminated when the interspike interval
exceeded 30 ms, Figure 1B. Thus defined, typical bursts had 5–17
spikes. Typical bursts lasted between 70 and 100 ms.
"""
spike_detect_burst_trn(spk; kwargs...) = spike_detect_burst(
    spk; t_silence=(0.070, 0.0), t_isi=0.030, nofs=(5, 0.07), kwargs...
)

@doc raw"""
    spike_detect_burst_lgn(spk; kwargs...)

Bursts were defined as two or more spikes, each spaced ≤4 ms apart following ≥100 ms of
silence, bursts rarely lasted more than 10 ms for LGN.
"""
spike_detect_burst_lgn(spk; kwargs...) = spike_detect_burst(
    spk; t_silence=(0.100, 0.0), t_isi=0.004, nofs=nothing, kwargs...
)

@doc raw"""
    spike_split_burst(spk; detector, kwargs...) -> (; burst::Vector{Vector}, tonic::Vector)

Split spike train into burst and tonic groups using `detector` function.
"""
function spike_split_burst(spk; detector=spike_detect_burst, kwargs...)
    _burst_idx = detector(spk; keep_index=true, kwargs...)
    _spk = deepcopy(spk)
    splice!(_spk, reduce(vcat, _burst_idx; init=Int64[]))
    return (; burst=map(x -> spk[x], _burst_idx), tonic=_spk)
end

@doc raw"""
    spike_split_burst_cardinal(spk; detector, kwargs...) -> (; burst::Vector, tonic::Vector)

Split spike train into burst and tonic groups using `detector` function.
But only the cardinal spike of burst is returned.
"""
function spike_split_burst_cardinal(spk; detector=spike_detect_burst, kwargs...)
    tmp = spike_split_burst(spk; detector, kwargs...)
    return (; tmp..., burst=map(first, tmp.burst))
end

@doc raw"""
    spike_burst_pattern(ISI::AbstractVector{<:AbstractVector{T}}; n=13, degree=Quadratic()) -> Matrix{T}

Interpolate burst inter-spike-intervals onto a common grid to generate comparable burst patterns.

## Arguments
- `ISI`: vector of ISI vectors, one per burst.
- `n`: number of interpolation points (default `13`).
- `degree`: B-spline degree for interpolation (default `Quadratic()`).

Returns an `n × length(ISI)` matrix.
"""
function spike_burst_pattern(ISI::AbstractVector{<:AbstractVector{T}}; n=13, degree=Quadratic()) where {T}
    output = zeros(T, n, length(ISI))
    Threads.@threads for i in eachindex(ISI)
        a = ISI[i]
        len = length(a)
        itp = interpolate(a, BSpline(degree))
        stride = (len - 1) / (n - 1)
        @inbounds for j in 1:n
            output[j, i] = itp(1 + (j - 1) * stride)
        end
    end
    return output
end
