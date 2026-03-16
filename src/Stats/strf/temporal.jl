## Temporal receptive fields
export trf_segment_biphasic, trf_segment_statistics
#TODO: export trf_isbiphasic

find_zero_crossing(trace) = findall(diff(signbit.(trace)) .!= 0)

function samesign_neighbor_until(trace, initnode, step=1)
    mysign = signbit(trace[initnode])
    
    N = length(trace)
    cursor = initnode + step
    
    while true
        ((cursor < 1) || (cursor > N)) && return cursor - step
    
        signbit(trace[cursor]) == mysign || return cursor - step
        cursor += step
    end
end

biphasic_rank_energy(trace, troi, proi) = sum(abs2, trace[[troi; proi]])
biphasic_rank_cor(trace, troi, proi) = begin
    tmp = zeros(size(trace))
    tmp[[troi; proi]] .= trace[[troi; proi]]
    cor(tmp, trace)
end

@doc raw"""
    trf_segment_biphasic(trf; by=biphasic_rank_energy)
This function will segment tSTA into a `triggering` phase and a `priming` phase.
TRF is assumed mean-subtracted already, that positive values will be treated as ON responses and negative values as OFF.
If `trf` does not have a `priming` phase, a `ArgumentError` will throw.

`by` option is setting how to compare different possible pairs of triggering/priming phases.
- `biphasic_rank_cor`: use overall correlation-coefficient with the raw trace (used by Ulas)
- `biphasic_rank_energy`: use squared sum (default, simpler and essentially same results)

It returns a NamedTuple of different ranges:
- triggering
- priming
"""
function trf_segment_biphasic(trf::AbstractVector; by=biphasic_rank_energy)
    # indices with zero crossing following
    # e.g., if it returns [2], then the crossing is between [2, 3].
    zero_crossing_ind0 = find_zero_crossing(trf) 

    ## no zero crossing found 
    isempty(zero_crossing_ind0) && throw(ArgumentError("no zero crossing found, please make sure it is z-scored or properly mean-shifted!"))

    ## found all putative triggering-priming pairs
    tp_pairs = []
    for ind0 in zero_crossing_ind0
        triggering_start = samesign_neighbor_until(trf, ind0, -1)
        triggering_end = ind0
        priming_start = ind0 + 1
        priming_end = samesign_neighbor_until(trf, priming_start, 1)

        push!(tp_pairs, (triggering_start:triggering_end, priming_start:priming_end))
    end

    ## pick one that has the largest ranking value
    (triggering_roi, priming_roi) = argmax(i -> by(trf, i...), tp_pairs)

    ## biological check
    last(triggering_roi) < 3 && @warn("triggering phase is too early ($(triggering_roi)), please check if it is biological.")

    return (; triggering=triggering_roi, priming=priming_roi)
end

@doc raw"""
    get_trf_statistics(trf, [segmentation]; kwargs...)

extract some basic statistics for each phases, including:
- `roi`
- `start` index
- `stop` index
- `peak` index
- `amplitude` value
- `magnitude` value
"""
function trf_segment_statistics(trf::AbstractVector, segmentation::NamedTuple; )
    N = length(trf)
    
    get_crossing_offset(trf, start, stop) = (stop-start) * trf[start] / (trf[start] - trf[stop])

    t0 = first(segmentation.triggering)
    t_start = t0 == 1 ? 1 : t0 + get_crossing_offset(trf, t0, t0 - 1)
    
    t1 = last(segmentation.triggering)
    p0 = first(segmentation.priming)
    zerocrossing = t1 + get_crossing_offset(trf, t1, p0)

    p1 = last(segmentation.priming)
    p_end = p1 == N ? N : p1 + get_crossing_offset(trf, p1, p1 + 1)
    
    (t_amplitude, t_peak) = findmax(abs, trf[segmentation.triggering])
    (p_amplitude, p_peak) = findmax(abs, trf[segmentation.priming])

    (;
        trace = trf,
        zerocrossing,
        triggering = (;
                    roi = segmentation.triggering,
                    start = t_start,
                    stop = zerocrossing,
                    peak = t0 + t_peak - 1,
                    amplitude = t_amplitude,
                    magnitude = sum(abs, trf[segmentation.triggering])
                    ),
        priming = (;
                    roi = segmentation.priming,
                    start = zerocrossing,
                    stop = p_end,
                    peak = p0 + p_peak - 1,
                    amplitude = p_amplitude,
                    magnitude = sum(abs, trf[segmentation.priming])
                    ),
    )
end

trf_segment_statistics(trf::AbstractVector; kwargs...) = trf_segment_statistics(trf, get_trf_segments_biphasic(trf; kwargs...))
