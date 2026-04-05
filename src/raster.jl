@doc raw"""
    spike_raster(spike_train::AbstractVector{T}, markers::AbstractVector; head=0.5, duration=1.0, tail=0.5, offset=true) -> Vector{Vector{T}}

Create rasters from `spike_train` and onset times of stimulus (`markers`),
with interval of `(marker - head, marker + duration + tail]`.

## Arguments
- `spike_train`: **sorted** vector of spike times.
- `markers`: onset times of each stimulus trial.
- `head`: pre-stimulus window (default `0.5` s).
- `duration`: stimulus duration (default `1.0` s).
- `tail`: post-stimulus window (default `0.5` s).
- `offset`: if `true` (default), spike times are relative to the onset time.
"""
function spike_raster(
    spk::AbstractVector,
    markers::AbstractVector;
    head::Real=0.5,
    duration::Real=1.0,
    tail::Real=0.5,
    offset::Bool=true,
)
    return map(markers) do item
        lo = searchsortedfirst(spk, item - head, lt=(<=))
        hi = searchsortedlast(spk, item + duration + tail)
        candidates = spk[lo:hi]
        offset ? candidates .- item : candidates
    end
end
