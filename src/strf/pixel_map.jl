@doc raw"""
    get_filtered_strf(strfs...; s=3)

Smooth each SRF with an `s × s` running-average kernel and normalise jointly by the
common absolute maximum.

## Arguments
- `strfs...`: one or more 2D spatial receptive field matrices.
- `s`: kernel size (default `3`).

Returns a tuple of smoothed, jointly normalised matrices.

See also: [`get_filtered_mask_full`](@ref), [`srf_pixel_cloud`](@ref).
"""
function get_filtered_strf(strfs...; s=3)
    strf_smoothed = map(i->conv(i, ones(s, s)./(s*s))[2:end-1, 2:end-1], strfs)
    strf_max = maximum(i -> maximum(abs, i), strf_smoothed)
    strf_smoothed_norm = map(i -> i ./ strf_max, strf_smoothed)
    return strf_smoothed_norm
end

@doc raw"""
    get_filtered_mask_full(strfs...; kwargs...)

Create binary masks of elements exceeding half the joint maximum.

## Arguments
- `strfs...`: one or more 2D SRF matrices.
- `kwargs...`: forwarded to `get_filtered_strf`.

In most cases, prefer [`srf_pixel_cloud`](@ref) which additionally applies flood-fill
connectivity.
"""
function get_filtered_mask_full(strfs... ; kwargs...)
    return map(i -> abs.(i) .> 0.5, get_filtered_strf(strfs...; kwargs...))
end

@doc raw"""
    srf_pixel_cloud(strfs::AbstractMatrix...; connectivity=4, s=3)

Filter and threshold SRFs at 50% of the joint maximum, then keep only
pixels connected to the peak pixel via flood fill.

## Arguments
- `strfs...`: 2D SRF matrices; processed collectively to find the joint maximum.
- `connectivity`: flood-fill connectivity, `4` or `8` (default `4`).
- `s`: moving-average kernel size (default `3`).

Returns a tuple of binary masks, one per input.

## Examples
```julia
(mask_on, mask_off) = srf_pixel_cloud(frame_on, frame_off)
(mask,) = srf_pixel_cloud(frame)
```

## References
- Müllner, F. E. & Roska, B. Individual thalamic inhibitory interneurons are functionally
  specialized toward distinct visual features. *Neuron* 112, 2765–2782.e9 (2024).
"""
function srf_pixel_cloud(strfs::AbstractMatrix...; connectivity=4, s=3)
    filtered_map = get_filtered_strf(strfs...; s)
    mask = get_filtered_mask_full(strfs...; s)

    initnode = let
        tmp = map(i->findmax(abs.(i)), filtered_map)
        argmax(i->i[1], tmp)[2]
    end

    _tmp = map(i->floodfill!(collect(Int8, i), initnode, 1, 9; connectivity), mask)
    return map(i->i .== 9, _tmp)
end
