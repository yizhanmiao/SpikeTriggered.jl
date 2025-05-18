export orientation_selectivity_index_global
export direction_selectivity_index_global

# function orientation_selectivity_index()
# end

# function direction_selectivity_index()
# end

@doc raw"""
    orientation_selectivity_index(resp::AbstractVector, theta::AbstractVector)

```math
    \text{gOSI} = \frac{\sum R_\theta e^{i 2\theta}}{\sum R_\theta}
```

### References:
- Gale, S. D. & Murphy, G. J. Distinct Representation and Distribution of Visual Information by Specific Cell Types in Mouse Superficial Superior Colliculus. J. Neurosci. 34, 13458-13471 (2014).
- Inayat, S. et al. Neurons in the Most Superficial Lamina of the Mouse Superior Colliculus Are Highly Selective for Stimulus Direction. J. Neurosci. 35, 7992-8003 (2015).
"""
function orientation_selectivity_index_global(resp::AbstractVector, theta::AbstractVector)
    _R = resp .* cis.(2 .* theta)
    return sum(_R) / sum(resp)
end

@doc raw"""
    direction_selectivity_index(resp::AbstractVector, theta::AbstractVector)

```math
    \text{gOSI} = \frac{\sum R_\theta e^{i \theta}}{\sum R_\theta}
```

### References:
- Gale, S. D. & Murphy, G. J. Distinct Representation and Distribution of Visual Information by Specific Cell Types in Mouse Superficial Superior Colliculus. J. Neurosci. 34, 13458-13471 (2014).
- Inayat, S. et al. Neurons in the Most Superficial Lamina of the Mouse Superior Colliculus Are Highly Selective for Stimulus Direction. J. Neurosci. 35, 7992-8003 (2015).
"""
function direction_selectivity_index_global(resp::AbstractVector, theta::AbstractVector)
    _R = resp .* cis.(theta)
    return sum(_R) / sum(resp)
end
