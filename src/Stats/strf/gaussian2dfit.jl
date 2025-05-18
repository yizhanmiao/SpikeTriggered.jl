export GaussianEllipse, make_ellipse
export srf_gaussian_fit

struct GaussianEllipse{T}
    amplitude::T
    axis_major::T
    axis_minor::T
    center_x::T
    center_y::T
    rotation::T
    bias::Union{T,Nothing}

    function GaussianEllipse(
        A::T, a::T, b::T, x::T, y::T, θ::T, bias::Union{T,Nothing}
    ) where {T}
        new{T}(A, a, b, x, y, θ, bias)
    end
    GaussianEllipse(param::AbstractArray{T}) where {T} = GaussianEllipse(param...)
end

area(e::GaussianEllipse; sigma=1) = abs2(sigma) * e.axis_major * e.axis_minor * π
function Base.collect(e::GaussianEllipse)
    _rez = [e.amplitude, e.axis_major, e.axis_minor, e.center_x, e.center_y, e.rotation]
    if !isnothing(e.bias)
        push!(_rez, e.bias)
    end
    return _rez
end

function make_ellipse(e::GaussianEllipse; kwargs...)
    _make_ellipse(e.axis_major, e.axis_minor, e.rotation, e.center_x, e.center_y; kwargs...)
end

function _make_ellipse(a, b, θ, x0, y0; σ=1, step=100)
    t = range(0; stop=2*pi, length=step)
    ellipse_x_r = a .* cos.(t)
    ellipse_y_r = b .* sin.(t)
    R = [cos(-θ) sin(-θ); -sin(-θ) cos(-θ)] .* σ
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = x0 .+ r_ellipse[:, 1]
    y = y0 .+ r_ellipse[:, 2]
    map(Point2, zip(x, y))
end

@doc raw"""
    image2dataset(image::AbstractMatrix) -> ([X Y], [Z])

convert 2d image into `(x, y)` coordinates and vector of elements (`z`).
Primarilly used for 2d gaussian fit.
"""
function image2dataset(image::AbstractMatrix{T}) where {T} # -> ([X Y], [Z])
    N = length(image)
    _coord = Matrix{Int64}(undef, N, 2)
    _z = Vector{T}(undef, N)
    for (idx, item) in enumerate(CartesianIndices(image))
        _coord[idx, 1] = item[2] # x
        _coord[idx, 2] = item[1] # y
        _z[idx] = image[item]
    end
    _coord, _z
end

function _loss_srf_gaussian_fit_least_square(X, Y)
    (param) -> begin
        (A, a, b, x0, y0, θ, bias) = length(param) == 6 ? [param; 0] : param

        â = cos(θ^2)/(2 * a^2) + sin(θ^2)/(2 * b^2)
        b̂ = sin(2*θ)/(4 * a^2) - sin(2*θ)/(4 * b^2)
        ĉ = sin(θ^2)/(2 * a^2) + cos(θ^2)/(2 * b^2)
        _deltaX = X[:, 1] .- x0
        _deltaY = X[:, 2] .- y0
        _raw = @. A * exp(-(â * _deltaX^2 + b̂ * _deltaX * _deltaY + ĉ * _deltaY^2)) +
            bias

        return mean(abs2, _raw .- Y)
    end
end

function srf_gaussian_fit(
    srf::AbstractMatrix{T}; param0=nothing, bias=true, optim_options=(;), raw=false
) where {T}
    X, Y = image2dataset(srf)
    param_0 = if isnothing(param0)
        _p = T[
            maximum(Y);
            1;
            1;
            mean(X, weights(Y), 1)[:];
            0.1;
        ]
        bias && push!(_p, mean(Y))
        _p
    else
        param0
    end

    loss = _loss_srf_gaussian_fit_least_square(X, Y)
    rez = optimize(loss, param_0, NelderMead(), OptimOptions(; optim_options...))
    return (raw ? rez : GaussianEllipse(rez.minimizer))
end

include("schiller_overlap.jl")
