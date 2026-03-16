export burst_detect, burst_detect_trn, burst_detect_lgn
export split_tonic_burst, split_tonic_cardinal
export burst_interpolate

@doc raw"""
    burst_detect(spk::Vector; t_silence=0.07, t_isi=0.03, nofs=nothing, keep_index=false) -> Vector{Vector{T}}

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
function burst_detect(
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
        output = Vector{T}[]
        output_index = Vector{Int64}[]
        for (burst, burst_index) in zip(burst_list, burst_index_list)
            if length(burst) >= nofs[1]
                if sum(diff(burst[1:nofs[1]])) <= nofs[2]
                    push!(output, burst)
                    push!(output_index, burst_index)
                end
            end
        end
        output, output_index
    else
        burst_list, burst_index_list
    end

    keep_index ? rez[2] : rez[1]
end

@doc raw"""
    detect_burst_trn(spk; kwargs...)

≥5 spikes within 70 ms, spaced ≤30 ms apart following ≥70 ms
of silence; bursts were terminated when the interspike interval
exceeded 30 ms, Figure 1B. Thus defined, typical bursts had 5–17
spikes. Typical bursts lasted between 70 and 100 ms.
"""
burst_detect_trn(spk; kwargs...) = burst_detect(
    spk; t_silence=(0.070, 0.0), t_isi=0.030, nofs=(5, 0.07), kwargs...
)

@doc raw"""
    detect_burst_lgn(spk; kwargs...)

Bursts were defined as two or more spikes, each spaced ≤4 ms apart following ≥100 ms of
silence, bursts rarely lasted more than 10 ms for LGN.
"""
burst_detect_lgn(spk; kwargs...) = burst_detect(
    spk; t_silence=(0.100, 0.0), t_isi=0.004, nofs=nothing, kwargs...
)

@doc raw"""
    split_tonic_burst(spk; detector, kwargs...) -> (; burst::Vector{Vector}, tonic::Vector)

Split spike train into burst and tonic groups using `detector` function.
"""
function split_tonic_burst(spk; detector=burst_detect, kwargs...)
    _burst_idx = detector(spk; keep_index=true, kwargs...)
    _spk = deepcopy(spk)
    splice!(_spk, reduce(vcat, _burst_idx; init=Int64[]))
    (; burst=map(x -> spk[x], _burst_idx), tonic=_spk)
end

@doc raw"""
    split_tonic_cardinal(spk; detector, kwargs...) -> (; burst::Vector, tonic::Vector)

Split spike train into burst and tonic groups using `detector` function.
But only the cardinal spike of burst is returned.
"""
function split_tonic_cardinal(spk; detector=burst_detect, kwargs...)
    tmp = split_tonic_burst(spk; detector, kwargs...)
    (; tmp..., burst=map(first, tmp.burst))
end

@doc raw"""
    interp_burst(burst_iti::Vector{Vector{T}; n=13, degree=:Quadratic()) -> Matrix{T} [n x nBursts]
    interp_burst(burst_iti::Vector{T}; n=13, interp_t=:steffen) -> Matrix{T} [n x 1]

Make interpolation of burst inter-spike-interval to generate burst patterns.

- `n` set the interpolation bin length.
- `degree` set the type of interpolation, which is using the `Interpolations`

The results will always be 2d matrix.
"""
function burst_interpolate(ISI::AbstractVector{<:AbstractVector{T}}; n=13, degree=Quadratic()) where {T}
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

@deprecate burst_interpolate(burst_iti::AbstractVector{<:Real}; kwargs...) burst_interpolate([burst_iti]; kwargs...)
