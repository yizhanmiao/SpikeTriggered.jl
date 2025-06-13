@doc raw"""
    spike_triggered_covariance(X::AbstractMatrix{T}, y::AbstractArray{T}; n::Integer=10, verbose=true) where {T, U} -> Matrix{T}
"""
function spike_triggered_covariance(
    X::AbstractMatrix{T}, y::AbstractArray{T}; n::Integer=10
) where {T}
    _Nx = size(X, 1)

    _sta = spike_triggered_average(X, y; n)
    _sta_c = _sta ./ sum(abs2, _sta)

    _stc = zeros(T, length(_sta), length(_sta))

    _y = mean(y, dims=2)[:]

    for idx in 1:length(_y)
        _ind = (((idx):(-1):(idx - n + 1)) .- 1) .% _Nx .+ 1
        _ind[_ind .< 1] .= _ind[_ind .< 1] .+ _Nx

        _s = X[_ind, :][:]
        _s = _s .- _s' * _sta .* _sta_c

        _stc .+= (_s * _s') .* y[idx]
    end
    _stc ./= sum(_y) - 1
    _stc
end
#TODO: review the function and make documentations
