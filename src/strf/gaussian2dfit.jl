@doc raw"""
    GaussianEllipse{T}

Ellipse representation of a 2D Gaussian fit result.

## Fields
- `amplitude`, `axis_major`, `axis_minor`, `center_x`, `center_y`, `rotation`
- `bias::Union{T,Nothing}` — optional additive bias term.

## Conversion
- `collect(e)` — parameter vector; `GaussianEllipse(vec)` reconstructs.

## Related functions
- [`area`](@ref), [`make_ellipse`](@ref), [`make_heatmap`](@ref),
  [`srf_schiller_overlap_index`](@ref).
"""
struct GaussianEllipse{T}
    amplitude::T
    axis_major::T
    axis_minor::T
    center_x::T
    center_y::T
    rotation::T
    bias::Union{T,Nothing}

    function GaussianEllipse(
        A::T, a::T, b::T, x::T, y::T, θ::T, bias::Union{T,Nothing}=nothing
    ) where {T}
        return new{T}(A, a, b, x, y, θ, bias)
    end
    GaussianEllipse(param::AbstractArray{T}) where {T} = GaussianEllipse(param...)
end

@doc raw"""
    area(e::GaussianEllipse; sigma=1) -> Real

Area of the ellipse at the given sigma contour.

```math
A = \sigma^2 \, \sigma_a \, \sigma_b \, \pi
```

## Arguments
- `e::GaussianEllipse`: ellipse to compute the area of.
- `sigma`: contour level in standard deviations (default `1`).
"""
area(e::GaussianEllipse; sigma=1) = abs2(sigma) * e.axis_major * e.axis_minor * π
function Base.collect(e::GaussianEllipse)
    _rez = [e.amplitude, e.axis_major, e.axis_minor, e.center_x, e.center_y, e.rotation]
    if !isnothing(e.bias)
        push!(_rez, e.bias)
    end
    return _rez
end

@doc raw"""
    make_ellipse(e::GaussianEllipse; reversed=true, σ=1, step=100) -> Vector{Point2}

Generate a polygon (vector of `Point2`) tracing the ellipse contour of a
[`GaussianEllipse`](@ref), suitable for plotting.

## Arguments
- `e::GaussianEllipse`: ellipse to trace.
- `reversed`: negate the rotation angle for image-coordinate systems where the
  y-axis points downward (default `true`).
- `σ`: contour level in standard deviations (default `1`).
- `step`: number of points on the polygon (default `100`).
"""
function make_ellipse(e::GaussianEllipse; reversed=true, kwargs...)
    #NOTE: rotation is inverted, because the y-axis is reversed during visualization
    θ = reversed ? e.rotation * -1 : e.rotation
    return _make_ellipse(e.axis_major, e.axis_minor, θ, e.center_x, e.center_y; kwargs...)
end

@doc raw"""
    make_heatmap(e::GaussianEllipse; xrange=1:16, yrange=nothing) -> Matrix

Evaluate the 2D Gaussian on a grid to produce a heatmap matrix.

## Arguments
- `e::GaussianEllipse`: Gaussian ellipse to evaluate.
- `xrange`: x-axis grid coordinates (default `1:16`).
- `yrange`: y-axis grid coordinates (default same as `xrange`).

Returns a matrix of size `length(yrange) × length(xrange)`.
"""
function make_heatmap(e::GaussianEllipse; xrange=1:16, yrange=nothing)
    (A, σx, σy, x0, y0, θ) = collect(e)[1:6]
    _func = gaussian_2d(; A, σx, σy, x0, y0, θ)
    yrange = isnothing(yrange) ? xrange : yrange
    return [_func(xi, yi) for yi in yrange, xi in xrange]
end

function _make_ellipse(a, b, θ, x0, y0; σ=1, step=100)
    t = range(0; stop=2 * pi, length=step)
    ellipse_x_r = a .* cos.(t)
    ellipse_y_r = b .* sin.(t)
    R = [cos(θ) sin(θ); -sin(θ) cos(θ)] .* σ
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = x0 .+ r_ellipse[:, 1]
    y = y0 .+ r_ellipse[:, 2]
    return map(Point2, zip(x, y))
end

@doc raw"""
    image2dataset(image::AbstractMatrix; upsample=1, interp=Linear(), norm=false) -> (X, Y)

Convert a 2D image into `(x, y)` coordinate–value pairs for 2D Gaussian fitting.

## Arguments
- `image`: 2D matrix to convert.
- `upsample`: upsampling ratio (default `1`, i.e. no upsampling).
- `interp`: interpolation method from `Interpolations.jl` (default `Linear()`).
- `norm`: normalise values by the absolute maximum (default `false`).

Returns `(X, Y)` where `X` is a matrix of `(x, y)` tuples and `Y` the corresponding values.
"""
function image2dataset(image::AbstractMatrix; upsample=1, interp=Linear(), norm=false)
    gridsize = size(image)
    _frm = interpolate(image, BSpline(interp))
    _sample_ind = [
        (xi, yi) for yi in 1:(1 / upsample):gridsize[2], xi in 1:(1 / upsample):gridsize[1]
    ]

    _x = map(identity, _sample_ind)
    _y = map(x -> _frm(x[2], x[1]), _sample_ind)
    #NOTE: # srf is organized as HxW, ie: [y, x]

    _y = norm ? _y ./ maximum(abs, _y) : _y

    return (_x, _y)
end

function gaussian_2d_param(; A, x0, y0, a, b, c, bias=0)
    return (x, y) -> begin
        dx = x - x0
        dy = y - y0
        dd = [dx; dy]

        rot = [
            a b
            b c
        ]
        A * exp(-sum(rot .* (dd * dd'))) + bias
    end
end

function gaussian_2d(; A, x0, y0, θ, σx, σy, bias=0)
    (a, b, c) = gaussian_2d_param_proj(; θ, σx, σy)
    return (x, y) -> begin
        A * exp(-(a * abs2(x - x0) + 2 * b * (x - x0) * (y - y0) + c * abs2(y - y0))) + bias
    end
end

function gaussian_2d_param_interp(; a, b, c)
    θ = atan(2b / (a - c)) / 2
    σx = 1 / sqrt(2(a * cos(θ)^2 + 2b * cos(θ)sin(θ) + c * sin(θ)^2))
    σy = 1 / sqrt(2(a * sin(θ)^2 - 2b * cos(θ)sin(θ) + c * cos(θ)^2))
    return (θ, σx, σy)
end

function gaussian_2d_param_proj(; θ, σx, σy)
    a = cos(θ)^2 / (2 * σx^2) + sin(θ)^2 / (2 * σy^2)
    b = -sin(θ)cos(θ) / (2 * σx^2) + sin(θ)cos(θ) / (2 * σy^2)
    c = sin(θ)^2 / (2 * σx^2) + cos(θ)^2 / (2 * σy^2)
    return (a, b, c)
end

function _loss_srf_gaussian_fit_least_square(X, Y)
    (param) -> begin
        (A, σx, σy, x0, y0, θ, bias) = length(param) == 6 ? [param; 0] : param
        _func = gaussian_2d(; A, σx, σy, x0, y0, θ, bias)
        ŷ = map(x -> _func(x...), X)
        return sum(abs2, ŷ .- Y)
    end
end

@doc raw"""
    srf_gaussian_fit(srf::AbstractMatrix{T}; kwargs...) -> GaussianEllipse{T}

Fit a 2D Gaussian to a spatial receptive field heatmap via least-squares optimisation.

## Arguments
- `srf`: 2D SRF matrix.
- `norm`: normalise by absolute maximum (default `false`).
- `upsample`: upsampling ratio (default `1`).
- `interp`: interpolation method from `Interpolations.jl` (default `Linear()`).
- `optim_algorithm`: `Optim.jl` algorithm (default `NelderMead()`).
- `optim_option`: `Optim.Options` for the solver.
- `param_init`: override default initial parameters.
- `bias`: include a bias term in the Gaussian (default `false`).
- `dev`: return the raw `Optim` result instead of a `GaussianEllipse` (default `false`).
- `verbose`: print optimisation output (default `false`).
"""
function srf_gaussian_fit(
    srf::AbstractMatrix{T};
    upsample=1,
    norm=false,
    interp=Linear(),
    optim_algorithm=NelderMead(),
    optim_option=OptimOptions(),
    param_init=nothing,
    bias=false,
    dev=false,
    verbose=false,
) where {T}
    (X, y) = image2dataset(srf; upsample, norm, interp)

    param_0 = if isnothing(param_init)
        _param = T[
            maximum(y), # A
            rand() * 2,
            rand() * 2, #σa, σb
            mean(first.(X), weights(y)), #x0
            mean(last.(X), weights(y)), # y0
            (rand() - 0.5) * pi / 2, # θ NOTE: init values within [-45deg, 45deg]
            # rand() # bias
        ]
        bias && push!(_param, rand(T)) # bias
        _param
    else
        collect(T, param_init)
    end

    _loss = _loss_srf_gaussian_fit_least_square(X, y)
    rez = optimize(_loss, param_0, optim_algorithm, optim_option)

    verbose && (@show rez)

    if dev
        return rez
    else
        _param = minimizer(rez)
        if length(_param) == 6
            return GaussianEllipse(_param..., nothing)
        else
            return GaussianEllipse(_param...)
        end
    end
end

include("schiller_overlap.jl")
# include("ellipse_overlap.jl")
