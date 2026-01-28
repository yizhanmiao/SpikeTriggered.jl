export get_pixel_cloud
public get_filtered_strf, get_filtered_mask_full

@doc raw"""
    get_filtered_strf(strfs...; s=3)

Perform a running average (3x3 by default, change `s` to have a different size)
on each strf input; and normalized by the common absolute maximum value.

### example
```julia
(s_on, s_off) = get_filtered_strf(frame_on, frame_off)
```

SEE ALSO: `get_filtered_mask_full` and `get_filtered_pixel_map`
"""
function get_filtered_strf(strfs...; s=3)
    strf_smoothed = map(i->conv(i, ones(s, s)./(s*s))[2:end-1, 2:end-1], strfs)
    strf_max = maximum(i -> maximum(abs, i), strf_smoothed)
    strf_smoothed_norm = map(i -> i ./ strf_max, strf_smoothed)
    return strf_smoothed_norm
end

@doc raw"""
    get_filtered_mask_full(strfs...)

return binary masks of elements larger than the half maximum.

In most cases, you should use `get_filtered_pixel_map`.
"""
function get_filtered_mask_full(strfs... ; kwargs...)
    return map(i -> abs.(i) .> 0.5, get_filtered_strf(strfs...; kwargs...))
end

@doc raw"""
    get_filtered_pixel_map(strfs...; connectivity=4, s=3)

Mimicing the behavior from [1], where on and off maps are filtered
and thresholded by 50% of the maximum; and only keep those pixels that
have direct connections with the maximum pixel.

The outputs are binary matrices.
                                         
## Arguments

- `strfs...`: frames of STAs; they will be used collectively to find the maximum value.
Make sure they are in 2D matrix format.

## Keyword Arguments

- `connectivity=4`: the final connection is using `flood fill` algorithm.
It can be either 4-way or 8-way connection.
- `s=3`: moving average kernel size.

## Examples

```julia
# process the on and off map together so the joint maximum value is used.
(mask_on, mask_off) = get_filtered_pixel_map(frame_on, frame_off)

# even with one input frame, the output is a one-element tuple.
(mask,) = get_filtered_pixel_map(frame)
```

1. Müllner, F. E. & Roska, B. Individual thalamic inhibitory interneurons are functionally specialized toward distinct visual features. Neuron 112, 2765-2782.e9 (2024).
"""
function get_pixel_cloud(strfs::AbstractMatrix...; connectivity=4, s=3)
    filtered_map = get_filtered_strf(strfs...; s)
    mask = get_filtered_mask_full(strfs...; s)
    
    initnode = let
        tmp = map(i->findmax(abs.(i)), filtered_map)
        argmax(i->i[1], tmp)[2]
    end

    _tmp = map(i->floodfill!(collect(Int8, i), initnode, 1, 9; connectivity), mask)
    return map(i->i .== 9, _tmp)
end

@deprecate get_filtered_pixel_map(args...; kwargs...) get_pixel_cloud(args...; kwargs...)
