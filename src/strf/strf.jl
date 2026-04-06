const AbstractSTRF{T} = AbstractArray{T,3} #[Height x Width x nTimepoints]
const AbstractSRF{T} = AbstractArray{T,2} #[Height x Width]
const AbstractTRF{T} = AbstractArray{T,1} #[nTimepoints]

@doc raw"""
    strf_(x, height, width) -> Array{T,3}
    strf_(x, gridsize=16) -> Array{T,3}
    strf_(x; gridsize=16) -> Array{T,3}

Reshape a flattened STA vector into a 3D spatiotemporal matrix (H × W × T).

## Arguments
- `x::AbstractArray`: flattened STA vector or matrix.
- `height`, `width`: spatial dimensions.
- `gridsize`: shorthand when height equals width (default `16`).

!!! note
    When plotting with `heatmap`, use the transposed matrix and reversed y-axis.
"""
strf_(x::AbstractArray, height::Integer, width::Integer) = reshape(x, height, width, :)
strf_(x::AbstractArray, gridsize::Integer) = strf_(x, gridsize, gridsize)
strf_(x::AbstractArray; gridsize=16) = strf_(x, gridsize, gridsize)

@doc raw"""
    strf_hstack(x::AbstractArray{T,3}) -> Matrix{T}
    strf_hstack(x, args...; kwargs...) -> Matrix{T}

Horizontally concatenate time slices of a 3D STRF into an H × (W × T) matrix.

## Arguments
- `x`: 3D STRF array, or arguments forwarded to `strf_` for reshaping first.

!!! note
    When plotting with `heatmap`, use the transposed matrix and reversed y-axis.
"""
strf_hstack(x::AbstractSTRF) = reduce(hcat, eachslice(x; dims=3))
strf_hstack(x, args...; kwargs...) = strf_hstack(strf_(x, args...; kwargs...))

@doc raw"""
    strf_vstack(strfs...; kwargs...) -> Matrix

Vertically stack multiple STRFs after horizontally concatenating each one.
So you can visualize multiple STRFs more compactively

## Arguments
- `strfs...`: one or more STRFs (3D arrays or arguments forwarded to `strf_hstack`).
- `kwargs...`: forwarded to `strf_hstack` (e.g. `gridsize`).
"""
strf_vstack(strfs...; kwargs...) = reduce(vcat, map(x->strf_hstack(x; kwargs...), strfs))

include("misc.jl")

# Spatial
include("gaussian2dfit.jl")
include("footprint.jl")
include("pixel_map.jl")

# Temporal
include("temporal.jl")
include("bright_dark_polarity.jl")
