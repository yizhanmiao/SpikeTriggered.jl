
@doc raw"""
    dispersion_index(psth; corrected=true, kwargs...) -> Real

Also known as the `Fano factor`.

```math
D = \frac{\sigma^2}{\mu}
```

## optional arguments:
- `dims`: calculate the index over dimension `dims`.

Check out more on [wikipedia](https://en.wikipedia.org/wiki/Index_of_dispersion).
"""
function dispersion_index(psth; corrected=true, kwargs...)
    return var(psth; corrected, kwargs...) ./ mean(psth; kwargs...)
end

# gaussian distribution cumulative density, single sided
# normcdf(x::T) where {T<:Real} = erfc(-x / sqrt(T(2))) / 2

# inverse of normcdf
# norminv(x::T) where {T<:Real} = erfcinv(2 * x) * (-sqrt(T(2)))

@doc raw"""
    modulation_index(Rpeak, Rnull) -> Real

Compute the modulation index (contrast ratio) between a peak and null response.

```math
MI = \frac{R_{\text{peak}} - R_{\text{null}}}{R_{\text{peak}} + R_{\text{null}}}
```

Returns a value in ``[-1, 1]``, where ``0`` indicates equal responses and
``\pm 1`` indicates complete dominance of one response.

## Arguments
- `Rpeak`: response at the preferred condition.
- `Rnull`: response at the null / opposite condition.
"""
function modulation_index(Rpeak, Rnull)
    return (Rpeak - Rnull) / (Rpeak + Rnull)
end

@doc raw"""
    troi2tproj(; head, duration, tail, step=0.002) -> StepRangeLen

Convert a time region of interest (TROI) into a uniformly spaced time-projection range.

## Arguments
- `head`: pre-stimulus time (positive value; the range starts at `-head`).
- `duration`: stimulus duration.
- `tail`: post-stimulus time.
- `step`: time step (default `0.002`, i.e. 2 ms).
"""
troi2tproj(; head, duration, tail, step=0.002) = range(-head; stop=duration + tail, step)
