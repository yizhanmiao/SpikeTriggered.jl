@doc raw"""
    osi_global(resp::AbstractVector, theta::AbstractVector)

Compute the global orientation selectivity index (gOSI).

```math
\text{gOSI} = \frac{\sum R_\theta e^{i 2\theta}}{\sum R_\theta}
```

## Arguments
- `resp`: response magnitude at each orientation.
- `theta`: orientation angles in radians.

Returns a complex number whose magnitude is the selectivity index.

## References
- Gale, S. D. & Murphy, G. J. Distinct Representation and Distribution of Visual
  Information by Specific Cell Types in Mouse Superficial Superior Colliculus.
  *J. Neurosci.* 34, 13458–13471 (2014).
- Inayat, S. et al. Neurons in the Most Superficial Lamina of the Mouse Superior
  Colliculus Are Highly Selective for Stimulus Direction. *J. Neurosci.* 35, 7992–8003 (2015).
"""
function osi_global(resp::AbstractVector, theta::AbstractVector)
    _R = resp .* cis.(2 .* theta)
    return sum(_R) / sum(resp)
end

@doc raw"""
    dsi_global(resp::AbstractVector, theta::AbstractVector)

Compute the global direction selectivity index (gDSI).

```math
\text{gDSI} = \frac{\sum R_\theta e^{i \theta}}{\sum R_\theta}
```

## Arguments
- `resp`: response magnitude at each direction.
- `theta`: direction angles in radians.

Returns a complex number whose magnitude is the selectivity index.

## References
- Gale, S. D. & Murphy, G. J. Distinct Representation and Distribution of Visual
  Information by Specific Cell Types in Mouse Superficial Superior Colliculus.
  *J. Neurosci.* 34, 13458–13471 (2014).
- Inayat, S. et al. Neurons in the Most Superficial Lamina of the Mouse Superior
  Colliculus Are Highly Selective for Stimulus Direction. *J. Neurosci.* 35, 7992–8003 (2015).
"""
function dsi_global(resp::AbstractVector, theta::AbstractVector)
    _R = resp .* cis.(theta)
    return sum(_R) / sum(resp)
end
