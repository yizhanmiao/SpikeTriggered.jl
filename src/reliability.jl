@doc raw"""
Schreiber reliability score:

```math
R_{\text{corr}} = \frac{2}{N(N-1)} \sum_{i = 1}^N \sum_{j = i+1} ^N \frac{\vec{s}_i \cdot \vec{s}_j}{|\vec{s}_i| |\vec{s}_j|}
```

Input in shape of `[nPSTH x nTrial]`

reference: dx.doi.org/10.1016/S0925-2312(02)00838-X
"""
function schreiber_score_matrix(psth::AbstractMatrix{T}) where {T<:Real}
    # [nPSTH x nTrial]
    N = size(psth, 2)
    _norm = sqrt.(sum(abs2, psth; dims=1))
    _norm = _norm' * _norm
    _dot = psth' * psth
    (sum(_dot ./ _norm) - N) / (N * (N - 1))
end

function schreiber_score_raster(
    raster::Vector{Vector{T}}; step=0.005, σ=0.003
) where {T<:Real}
    _start = floor(minimum(x->minimum(x; init=Inf), raster) / step)
    _end = ceil(maximum(x->maximum(x; init=(-Inf)), raster) / step)
    _trange = collect(T, range(_start; stop=_end, step))
    _psth = map(x->spike_histogram_smoothed(x, _trange; σ), raster)
    schreiber_score_matrix(reduce(hcat, _psth))
end
