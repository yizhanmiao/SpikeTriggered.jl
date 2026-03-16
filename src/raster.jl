export spike_raster

@doc raw"""
    spike_raster(spike_train::AbstractVector{T}, markers::AbstractVector; head=0.5, duration=1.0, tail=0.5, offset=true) -> Vector{Vector{T}}

Create rasters from `spike_train` and onset times of stimulus (`markers`),
with interval of `(marker-head, marker+duration+tail]`.

If `offset` is `true`, spike times of each trial be relative to the onset time.
"""
function spike_raster(
    spk::AbstractVector{T},
    markers::AbstractVector;
    head::Real=0.5,
    duration::Real=1.0,
    tail::Real=0.5,
    offset::Bool=true,
) where {T}
    output = Vector{T}[]
    for item in markers
        _candidates = spk[(item - head) .< spk .<= (item + duration + tail)]
        if offset
            push!(output, _candidates .- item)
        else
            push!(output, _candidates)
        end
    end
    return output
end
