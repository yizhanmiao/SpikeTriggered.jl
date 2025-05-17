export srf_schiller_overlap_index

@doc raw"""
    srf_schiller_overlap_index(e1::GaussianEllipse, e2::GaussianEllipse)

overlap of two spatial receptive fields based on Schiller's overlap index.

For two receptive fields with the Gaussian parameters ($\alpha$, $\sigma_a$, $\sigma_b$, $\theta$, $x_0$, $y_0$), the 2D overlap index is defined as:

```math
\Omega = \frac{\sigma_{12} + \sigma_{21} - d}{\sigma_{12} + \sigma_{21} + d}
```

Where $d$ is the Euclidian distance between the centers, and:

```math
\sigma_{ij} = \frac{\sigma_{ai}\sigma_{bi}}{\sqrt{
\sigma_{ai}^2 \sin^2(\theta_i - \phi_{ij}) +
\sigma_{bi}^2 \cos^2(\theta_i - \phi_{ij})
}}
```

```math
\phi_{ij} = \arctan\left(\frac{y_{0j} - y_{0i}}{x_{0j} - x_{0i}}\right)
```

### References
- Schiller et al., 1979
- Wang et al., 2007
"""
function srf_schiller_overlap_index(e1::GaussianEllipse, e2::GaussianEllipse)
    _sigma_sum = _schiller_sigma(e1, e2) + _schiller_sigma(e2, e1)
    _d = sum(abs2, (e1.center_x - e2.center_x, e1.center_y - e2.center_y)) |> sqrt
    (_sigma_sum - _d) / (_sigma_sum + _d)
end

function _schiller_sigma(ei, ej)
    _phi = _schiller_phi(ei, ej)
    ei.axis_major * ei.axis_minor / sqrt(
        (ei.axis_major * sin(ei.rotation - _phi))^2 +
        (ei.axis_minor * cos(ei.rotation - _phi))^2
    )
end

function _schiller_phi(ei, ej)
    atan((ej.center_y-ei.center_y)/(ej.center_x-ei.center_x))
end