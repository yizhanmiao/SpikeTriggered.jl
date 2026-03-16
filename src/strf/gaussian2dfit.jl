export GaussianEllipse, area, make_ellipse, make_heatmap
export srf_gaussian_fit

@doc raw"""
    GaussianEllipse{T}

Ellipse representation for 2D gaussian fit result

## Conversion:
- `collect(e::GaussianEllipse)`: convert ellipse to a vector
- `GaussianEllipse(param::Vector)`: convert the vector above back to an ellipse.

## Functions:
- `area(e::GaussianEllipse; sigma=1)`: get the ellipse area for a given sigma value.
- `make_ellipse(e::GaussianEllipse)`: get the contour line for plotting
- `make_heatmap(e::GaussianEllipse)`: get the 2d heatmap for plotting
- `srf_schiller_overlap_index(e1, e2)`
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

area(e::GaussianEllipse; sigma=1) = abs2(sigma) * e.axis_major * e.axis_minor * π
function Base.collect(e::GaussianEllipse)
    _rez = [e.amplitude, e.axis_major, e.axis_minor, e.center_x, e.center_y, e.rotation]
    if !isnothing(e.bias)
        push!(_rez, e.bias)
    end
    return _rez
end

function make_ellipse(e::GaussianEllipse; reversed=true, kwargs...)
    #NOTE: rotation is inverted, because the y-axis is reversed during visualization
    θ = reversed ? e.rotation * -1 : e.rotation
    return _make_ellipse(e.axis_major, e.axis_minor, θ, e.center_x, e.center_y; kwargs...)
end

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
    image2dataset(image::AbstractMatrix; upsample=1, interp=Linear(), norm=false) -> ([X Y], [Z])

convert 2d image into matrix of `(x, y)` coordinates and matrix of elements (`z`).
Primarilly used for 2d gaussian fit.

Keyword Arguments:
- upsample: ratio of upsampling, default as `1`` (i.e. no upsampling)
- interp: upsampling method, default as `Linear()`` (using the `Interpolations.jl` package)
- norm: flag to normalize the image value by its maximum, default as `false`
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
    (θ, σx, σy)
end

function gaussian_2d_param_proj(; θ, σx, σy)
    a = cos(θ)^2 / (2 * σx^2) + sin(θ)^2 / (2 * σy^2)
    b = -sin(θ)cos(θ) / (2 * σx^2) + sin(θ)cos(θ) / (2 * σy^2)
    c = sin(θ)^2 / (2 * σx^2) + cos(θ)^2 / (2 * σy^2)
    (a, b, c)
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
    srf_gaussian_fit(
        srf::AbstractMatrix{T};
        norm=false, upsample=1, interp=Linear(),
        optim_algorithm=NelderMead(),
        optim_option=OptimOptions(),
        param_init=nothing,
        bias=false,
        dev=false,
        verbose=false,
    ) -> GaussianEllipse{T}

2D Gaussian fit of one spatial receptive field heatmap.

Keyword Arguments:
- norm: flag to normalize the image value by its maximum, default as `false`
- upsample: ratio of upsampling, default as `1`` (i.e. no upsampling)
- interp: upsampling method, default as `Linear()`` (using the `Interpolations.jl` package)
- optim_algorithm: optimization algorithm from `Optim.jl`, default as `NelderMead()` (which should be the same algorithm of `fminsearch` in MATLAB)
- optim_option: Optim.Options
- param_init: initial parameters to overwrite the default method.
- bias: flag to include a bias term in the 2D Gaussian function, default as `false`.
- dev: flag to return the Optim result, otherwise return a GaussianEllipse object, default as `false`.
- verbose: to print optimization output, default as `false`.
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
        rez
    else
        _param = minimizer(rez)
        if length(_param) == 6
            GaussianEllipse(_param..., nothing)
        else
            GaussianEllipse(_param...)
        end
    end
end

include("schiller_overlap.jl")
# include("ellipse_overlap.jl")
