@doc raw"""
    spike_train_spectrum_power(spks::AbstractSpikeTrain{T}, freq; t_len=1) where {T <: AbstractFloat}

The spike-train power spectrum is computed from the Fourier transform of the spike train,
as in (Dummer, et al., 2014).

```math
S(f) = \frac{\langle \tilde{x} \tilde{x}^\ast \rangle}{T}
```

The neural spike train is represented by a sum of delta functions at the spike times:
```math
x(t) = \sum_i \delta(t - t_i)
```

The Fourier transform for the time window is defined by
```math
\tilde{x}(f) = \int_0^{T} e^{- 2 \pi j f t} x(t)\; dt
```

### Reference:
- Dummer, B., Wieland, S. & Lindner, B. Self-consistent determination of the spike-train power spectrum in a neural network with sparse connectivity. Front. Comput. Neurosci. 8, (2014). DOI: 10.3389/fncom.2014.00104

"""
function spike_train_spectrum_power(spks::AbstractSpikeTrain{T}, f; t_len=1) where {T <: AbstractFloat}
    ω = T(2 * pi * f)
    ϕ = spks .* ω
    _real = sum(cos, ϕ)
    _imag = sum(sin, ϕ)
    return (_real^2 + _imag^2) / t_len
end