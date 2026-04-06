## Bright-Dark Polarity Score
# it is represented as a complex number,
# Omega = Ron + Roff * im = R * cis(theta)
#
# R = abs(Omega)
# theta = angle(Omega)
# Ron = real(Omega)
# Roff = imag(Omega)

@doc raw"""
    trf_polarity_score(data; thresh=0.8, reverse_off=true)
    trf_polarity_score(on, off; kwargs...)

Compute the bright–dark polarity score, a complex-valued index that encodes the
relative strength of ON and OFF responses at the triggering latency.

```math
\Omega = R_\text{on} + R_\text{off} \, i = r \, e^{i\theta}
```

The angle ``\theta`` classifies the cell type:
- ``(-\pi,\, -\pi/2)``: `:sbc` (suppressed-by-contrast)
- ``[-\pi/2,\, 0]``: `:on`
- ``(0,\, \pi/2)``: `:onoff`
- ``[\pi/2,\, \pi]``: `:off`

## Arguments
- `data`: a NamedTuple (or similar) with fields `:on` and `:off`, each a TRF vector.
  Alternatively, pass `on` and `off` vectors directly.
- `thresh`: z-score threshold for detecting the first non-zero derivative (default `0.8`).
- `reverse_off`: negate the OFF response before combining (default `true`).

Returns `(; triggering_latency, bdscore, bdtheta, bdtype)`.
"""
function trf_polarity_score(data; thresh=0.8, reverse_off=true)
    if !(haskey(data, :on) && haskey(data, :off))
        throw(ArgumentError("data must be a NamedTuple with fields `on` and `off`"))
    end
    trf = (data.on, data.off)

    trf_der = map(trf) do item
        _der = diff(item)
        _der_z = [0; (_der .- mean(_der)) ./ std(_der)]
        _der_z[abs.(_der_z) .< thresh] .= 0
        _der_z
    end

    triggering_latency = minimum(trf_der) do item
        ind = findfirst(!iszero, item)
        isnothing(ind) ? typemax(Int64) : ind
    end

    bdscore = if triggering_latency == typemax(Int64)
        0 + 0 * im
    else
        Ron = data.on[triggering_latency]
        Roff = reverse_off ? -1 * data.off[triggering_latency] : data.off[triggering_latency]
        Ron + Roff * im
    end

    bdtheta = angle(bdscore)
    bdtype = if -pi < bdtheta < -pi/2
        :sbc
    elseif -pi/2 <= bdtheta <= 0
        :on
    elseif 0 < bdtheta < pi/2
        :onoff
    elseif pi/2 <= bdtheta <= pi
        :off
    else
        :unknown
    end

    return (;
        triggering_latency,
        bdscore,
        bdtheta=rad2deg(bdtheta),
        bdtype,
    )
end

trf_polarity_score(on::AbstractTRF, off::AbstractTRF; kwargs...) = trf_polarity_score((; on, off); kwargs...)
