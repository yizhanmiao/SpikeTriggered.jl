export strf_, strf_hstack, strf_vstack

const STRF{T} = Array{T,3} #[Height x Width x nTimepoints]
const SRF{T} = Array{T,2} #[Height x Width]
const TRF{T} = Array{T,1} #[nTimepoints]

@doc raw"""
    strf_(val, height, width) -> Array{T, 3} HxWxT
    strf_(val, gridsize=16) -> Array{T, 3} HxWxT
    strf_(val; gridsize=16) -> Array{T, 3} HxWxT

convert strf vector into 3d spatiotemporal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
strf_(x::AbstractArray, height::Integer, width::Integer) = reshape(x, height, width, :)
strf_(x::AbstractArray, gridsize::Integer) = strf_(x, gridsize, gridsize)
strf_(x::AbstractArray; gridsize=16) = strf_(x, gridsize, gridsize)

@doc raw"""
    strf_hstack(strf, args...; kwargs...) -> Array{T, 2} Sx(SxT)

convert strf 3d matrix into horizontal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
strf_hstack(x::AbstractArray{T,3}) where {T} = reduce(hcat, eachslice(x; dims=3))
strf_hstack(x, args...; kwargs...) = strf_hstack(strf_(x, args...; kwargs...))

@doc raw"""
    strf_vstack(strfs...; kwargs...)

stack multiple strfs vertically.
"""
strf_vstack(strfs...; kwargs...) = reduce(vcat, map(x->strf_hstack(x; kwargs...), strfs))

@deprecate vstack_strf strf_hstack
@deprecate hstack_strf strf_hstack
@deprecate make_strf strf_

include("gaussian2dfit.jl")
include("misc.jl")
include("pixel_map.jl")
