@doc raw"""
    gaussian_filter_gsl(arr::Vector{Float64}; K::Integer, α::Real=1)

The GSL gaussian kernel is defined by:
```math
G(k) = \exp\left(- \frac{k^2}{2 \sigma^2}\right)
```

for ``-(K-1)/2 \leq k \leq (K-1)/2``, and ``K`` is the size of the kernel.
The parameter ``\alpha`` specifies the number of standard deviations ``\sigma``
desired in the kernel.
```math
\sigma = \frac{K-1}{2 \alpha}
```
"""
function gaussian_filter_gsl(arr::Vector{Float64}; K::Integer, α::Real=1)
    _output = zeros(length(arr))

    _x = GSL.C.gsl_vector(length(arr), 1, pointer(arr), C_NULL, 0)
    _y = GSL.C.gsl_vector(length(_output), 1, pointer(_output), C_NULL, 0)

    _workspace = GSL.filter_gaussian_alloc(K)
    GSL.filter_gaussian(GSL.GSL_FILTER_END_PADZERO, α, 0, _x, _y, _workspace)
    GSL.filter_gaussian_free(_workspace)

    return _output
end

function _gaussian_kernel(::Type{T}, K::Integer, α::Real; norm=false) where {T}
    output = zeros(K)
    N = div(K - 1, 2)
    _x = (-N):N
    output = T.(exp.(-0.5 .* (α .* _x ./ N) .^ 2))
    if norm
        output ./ sum(output)
    else
        output
    end
end

function filter_gaussian(
    arr::Vector{T}; K, α, norm=true, norm_tail=false
) where {T<:AbstractFloat}
    _kernel = _gaussian_kernel(T, K, α; norm)
    _n = div(K - 1, 2)
    _filted = conv(arr, _kernel)[(_n + 1):(end - _n)]

    output = if norm_tail
        _track = conv(ones(T, length(arr)), _kernel)[(_n + 1):(end - _n)]
        _filted ./ _track
    else
        _filted
    end
end

function filter_LP_butterworth(trace::Vector{T}; cutoff, order=4, fs=10000) where {T}
    responsetype = Lowpass(cutoff; fs)
    designmethod = Butterworth(order)

    filtfilt(digitalfilter(responsetype, designmethod), trace)
end

function filter_BP_butterworth(trace::Vector{T}; cutoff, order=4, fs=10000) where {T}
    responsetype = Bandpass(cutoff...; fs)
    designmethod = Butterworth(order)

    filtfilt(digitalfilter(responsetype, designmethod), trace)
end

function filter_HP_butterworth(trace::Vector{T}; cutoff, order=4, fs=10000) where {T}
    responsetype = Highpass(cutoff; fs)
    designmethod = Butterworth(order)

    filtfilt(digitalfilter(responsetype, designmethod), trace)
end

function filter_median(x::Vector{T}; order=9) where {T<:Integer}
    @assert(order % 2 == 1, "order should be odd number.")

    N = length(x)
    output = zeros(T, N)
    _width = div(order - 1, 2)

    for idx in 1:_width
        @inbounds output[idx] = round(T, median(x[1:idx]))
    end

    Threads.@threads for idx in (_width + 1):(N - _width)
        @inbounds output[idx] = median(x[(idx - _width):(idx + _width)])
    end

    for idx in (N - _width + 1):N
        @inbounds output[idx] = round(T, median(x[idx:N]))
    end

    output
end

function filter_median(x::Vector{T}; order=9) where {T<:Real}
    @assert(order % 2 == 1, "order should be odd number.")

    N = length(x)
    output = zeros(T, N)
    _width = div(order - 1, 2)

    for idx in 1:_width
        @inbounds output[idx] = median(x[1:idx])
    end

    Threads.@threads for idx in (_width + 1):(N - _width)
        @inbounds output[idx] = median(x[(idx - _width):(idx + _width)])
    end

    for idx in (N - _width + 1):N
        @inbounds output[idx] = median(x[idx:N])
    end

    output
end
