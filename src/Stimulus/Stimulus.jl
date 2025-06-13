module Stimulus

using Base: Base

include("sparse_noise.jl")
include("markers.jl")

@doc raw"""
    apply_video_nonlinearity(video, bias=0) -> video

For sparse noise stimuli, to retrieve nonlinear map (ON map and OFF map)
by ignore bright/dark spots.

# Arguments
- video: the original video matrix
- bias: could be either values or symbols
    - ON: `:on` or `-1`
    - OFF: `:off` or `1`
    - DIFF: `:diff` or `0`
"""
function apply_video_nonlinearity(video, bias=0)
    _bias = if isa(bias, Symbol)
        if bias == :on
            -1
        elseif bias == :off
            1
        else # e.g. :diff
            0
        end
    else
        bias
    end
    _video = copy(video)
    _video[_video .== _bias] .= 0
    _video
end

end
