export event_triggered_raster

function event_triggered_raster(spike, marker, marker_mask; kwargs...)
    _raster = spike_raster(spike, marker; kwargs...)
    return event_triggered_raster(_raster, marker_mask)
end

function event_triggered_raster(
    raster::SpikeRaster{U}, marker_mask::AbstractVector{T};
) where {T,U}
    @assert length(raster) == length(marker_mask) "length of raster and marker_mask must be equal."
    _unique_ids = sort(unique(marker_mask))
    output = Vector{Pair{T,SpikeRaster{U}}}(undef, length(_unique_ids))
    for (idx, i) in enumerate(_unique_ids)
        _roi = marker_mask .== i
        output[idx] = i => raster[_roi]
    end
    return output
end

#TODO: review the function and make documentations
