# Statistics

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
    var(psth; corrected, kwargs...) ./ mean(psth; kwargs...)
end

# gaussian distribution cumulative density, single sided
normcdf(x::T) where {T<:Real} = erfc(- x / sqrt(T(2))) / 2

# inverse of normcdf
norminv(x::T) where {T<:Real} = erfcinv(2 * x) * (- sqrt(T(2)))

function modulation_index(Rpeak, Rnull)
    (Rpeak - Rnull) / (Rpeak + Rnull)
end
