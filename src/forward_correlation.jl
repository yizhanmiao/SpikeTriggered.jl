

@doc raw"""
    spike_raster_groupby(spike_train, marker, marker_mask; kwargs...) -> Vector{Pair{T, SpikeRaster}}

Compute a spike raster and regroup trials by stimulus condition.

Calls `spike_raster` with the given parameters, then groups the resulting
trials by their `marker_mask` labels.

## Arguments
- `spike_train::AbstractSpikeTrain`: **sorted** vector of spike times.
- `marker::AbstractVector`: stimulus onset times.
- `marker_mask::AbstractVector`: condition label for each trial (integers recommended).
- `kwargs...`: forwarded to `spike_raster` (e.g. `head`, `duration`, `tail`, `offset`).

Returns a sorted vector of `label => SpikeRaster` pairs.
"""
function spike_raster_groupby(
    spike_train::AbstractSpikeTrain,
    marker::AbstractVector,
    marker_mask::AbstractVector;
    kwargs...,
)
    _raster = spike_raster(spike_train, marker; kwargs...)
    return spike_raster_groupby(_raster, marker_mask)
end

@doc raw"""
    spike_raster_groupby(raster, marker_mask) -> Vector{Pair{T, SpikeRaster}}

Regroup a pre-computed raster into distinct stimulus-condition groups.
This is essentially a forward correaltion.

## Arguments
- `raster::SpikeRaster`: collection of spike trains across trials.
- `marker_mask::AbstractVector`: condition label for each trial; must have the same
  length as `raster`. Integers are recommended for clear grouping, but it can be anything.

Returns a sorted vector of `label => SpikeRaster` pairs, one entry per unique label.
"""
function spike_raster_groupby(
    raster::SpikeRaster{U}, marker_mask::AbstractVector{T}
) where {T,U}
    if length(raster) != length(marker_mask)
        throw(ArgumentError("raster length ($(length(raster))) must equal marker_mask length ($(length(marker_mask)))"))
    end
    _unique_ids = sort(unique(marker_mask))
    return map(_unique_ids) do i
        i => raster[marker_mask .== i]
    end
end
