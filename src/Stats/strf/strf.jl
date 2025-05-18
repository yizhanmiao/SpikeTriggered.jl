export strf_, strf_hstack

const STRF{T} = Array{T,3} #[spatialX x spatialY x nTimepoints]
const SRF{T} = Array{T,2} #[spatialX x spatialY]
const TRF{T} = Array{T,1} #[nTimepoints]

strf_(x::AbstractArray, width::Integer, height::Integer) = reshape(x, width, height, :)
strf_(x::AbstractArray; gridsize=16) = strf_(x, gridsize, gridsize)
strf_(x::AbstractArray{T,2}) where {T} = begin
    nd = size(x, 1)
    gridsize = sqrt(nd)
    if isinteger(gridsize)
        strf_(x, Int(gridsize), Int(gridsize))
    else
        error("Input array is not a square matrix")
    end
end

strf_hstack(x::AbstractArray{T,3}) where {T} = reduce(hcat, eachslice(x; dims=3))
strf_hstack(x, args...; kwargs...) = strf_hstack(strf_(x, args...; kwargs...))

@doc raw"""
    strf_(val, width, height) -> Array{T, 3} XxYxT
    strf_(val; gridsize=16) -> Array{T, 3} XxYxT
    strf_(val::Matrix) -> Array{T, 3} XxYxT

convert strf vector into 3d spatiotemporal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
strf_

@doc raw"""
    strf_hstack(strf) -> Array{T, 2} Sx(SxT)
    strf_hstack(strf; gridsize) -> Array{T, 2} Sx(SxT)
    strf_hstack(strf, width, height) -> Array{T, 2} Sx(SxT)

convert strf 3d matrix into horizontal matrix.

When plotting with heatmap, REMEMBER to use
the transversed matrix and reversed y axis.
"""
strf_hstack

@deprecate hstack_strf strf_hstack
@deprecate make_strf strf_

include("gaussian2dfit.jl")
