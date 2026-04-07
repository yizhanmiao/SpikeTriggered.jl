## Temporal receptive fields

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
    trf_segment_biphasic(trf::AbstractVector; by=biphasic_rank_energy)

Segment a temporal receptive field into a `triggering` phase and a `priming` phase.
The TRF is assumed to be mean-subtracted (positive = ON, negative = OFF).

## Arguments
- `trf`: mean-subtracted temporal receptive field vector.
- `by`: ranking function to select the best triggering/priming pair.
  - `biphasic_rank_energy` — squared sum (default).
  - `biphasic_rank_cor` — correlation with the raw trace.

Returns `(; triggering, priming)` where each value is an index range.

Throws `ArgumentError` if no zero crossing is found.
"""
function trf_segment_biphasic(trf::AbstractVector; by=biphasic_rank_energy)
    # indices with zero crossing following
    # e.g., if it returns [2], then the crossing is between [2, 3].
    zero_crossing_ind0 = find_zero_crossing(trf)

    ## no zero crossing found
    isempty(zero_crossing_ind0) && throw(
        ArgumentError(
            "no zero crossing found, please make sure it is z-scored or properly mean-shifted!",
        ),
    )

    ## found all putative triggering-priming pairs
    tp_pairs = map(zero_crossing_ind0) do ind0
        triggering_start = samesign_neighbor_until(trf, ind0, -1)
        priming_start = ind0 + 1
        priming_end = samesign_neighbor_until(trf, priming_start, 1)
        (triggering_start:ind0, priming_start:priming_end)
    end

    ## pick one that has the largest ranking value
    (triggering_roi, priming_roi) = argmax(i -> by(trf, i...), tp_pairs)

    ## biological check
    last(triggering_roi) < 3 && @warn(
        "triggering phase is too early ($(triggering_roi)), please check if it is biological."
    )

    return (; triggering=triggering_roi, priming=priming_roi)
end

@doc raw"""
    trf_segment_statistics(trf::AbstractVector, segmentation::NamedTuple)
    trf_segment_statistics(trf::AbstractVector; kwargs...)

Extract statistics for each phase of a biphasic TRF segmentation.

## Arguments
- `trf`: temporal receptive field vector.
- `segmentation`: output of `trf_segment_biphasic` (or computed automatically via `kwargs`).

Returns a named tuple with `trace`, `zerocrossing`, and per-phase statistics
(`roi`, `start`, `stop`, `peak`, `amplitude`, `magnitude`) for `triggering`
and `priming`.
"""
function trf_segment_statistics(trf::AbstractVector, segmentation::NamedTuple;)
    N = length(trf)

    get_crossing_offset(trf, start, stop) =
        (stop-start) * trf[start] / (trf[start] - trf[stop])

    t0 = first(segmentation.triggering)
    t_start = t0 == 1 ? 1 : t0 + get_crossing_offset(trf, t0, t0 - 1)

    t1 = last(segmentation.triggering)
    p0 = first(segmentation.priming)
    zerocrossing = t1 + get_crossing_offset(trf, t1, p0)

    p1 = last(segmentation.priming)
    p_end = p1 == N ? N : p1 + get_crossing_offset(trf, p1, p1 + 1)

    (t_amplitude, t_peak) = findmax(abs, trf[segmentation.triggering])
    (p_amplitude, p_peak) = findmax(abs, trf[segmentation.priming])

    return (;
        trace=trf,
        zerocrossing,
        triggering=(;
            roi=segmentation.triggering,
            start=t_start,
            stop=zerocrossing,
            peak=t0 + t_peak - 1,
            amplitude=t_amplitude,
            magnitude=sum(abs, trf[segmentation.triggering]),
        ),
        priming=(;
            roi=segmentation.priming,
            start=zerocrossing,
            stop=p_end,
            peak=p0 + p_peak - 1,
            amplitude=p_amplitude,
            magnitude=sum(abs, trf[segmentation.priming]),
        ),
    )
end

function trf_segment_statistics(trf::AbstractVector; kwargs...)
    trf_segment_statistics(trf, trf_segment_biphasic(trf; kwargs...))
end
