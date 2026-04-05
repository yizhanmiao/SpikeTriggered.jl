@doc raw"""
    srf_footprint([func_pval,] raw::AbstractArray; alpha=0.01, fdr_c=:auto, collapse=true, gridsize=nothing)

Create a signed footprint map via Benjamini–Hochberg FDR control on q-values
from `func_pval(raw)`.

```math
P_k \leq \frac{k}{N \cdot C} \alpha
```

For arbitrary joint distributions of p-values, `fdr_c` defaults to
`benjamini_hochberg_constant(N)` (≈ ln(N) + γ + 1/N).

## Arguments
- `func_pval`: function mapping raw values to p-values (default `zscore_pvalue`).
- `raw`: array of test statistics.
- `alpha`: FDR threshold (default `0.01`).
- `fdr_c`: BH constant (default `:auto`).
- `collapse`: return a 2D time-collapsed map (default `true`).
- `gridsize`: spatial grid size for collapsing.

Sign codes: `1` ON, `0` no response, `-1` OFF.
"""
function srf_footprint(
    func_pval::Function, raw::AbstractArray; collapse=true, gridsize=nothing, kwargs...
)
    gridsize = if collapse && gridsize isa Nothing
        if ndims(raw) > 1
            size(raw, 1)
        else
            throw(ArgumentError("please specify gridsize to get collapsed footprint map"))
        end
    else
        gridsize
    end

    _mask = srf_footprint_mask(func_pval, raw; kwargs...)
    _map = sign.(_mask .* raw)

    if collapse
        return collapse_footprint(_map; gridsize)
    else
        return _map
    end
end

function srf_footprint(raw::AbstractArray; kwargs...)
    return srf_footprint(zscore_pvalue, raw; kwargs...)
end

@doc raw"""
    collapse_footprint(fp; kwargs...)

Collapse a spatiotemporal footprint into a spatial footprint. Polarity is
determined by the first non-zero sign along the temporal axis.

## Arguments
- `fp`: flattened footprint array.
- `kwargs...`: forwarded to `strf_` for reshaping (e.g. `gridsize`).
"""
function collapse_footprint(fp; kwargs...)
    return map(
        #NOTE: polarity is set as the first sign of footprint
        x -> begin
            _idx = findfirst(abs.(x) .> 0)
            _idx isa Nothing ? 0 : x[_idx]
        end,
        eachslice(strf_(fp; kwargs...); dims=(1, 2)),
    )
end

@doc raw"""
    srf_footprint_mask([func_pval::Function,] raw::AbstractArray; alpha=0.01, fdr_c=:auto)

Create a binary footprint mask via Benjamini–Hochberg FDR control.

## Arguments
- `func_pval`: function mapping raw values to p-values (default `zscore_pvalue`,
  which assumes standard-normal z-scores).
- `raw`: array of test statistics.
- `alpha`: FDR threshold (default `0.01`).
- `fdr_c`: BH constant for correcting joint p-value distributions (default `:auto`).

In most cases, prefer [`srf_footprint`](@ref) which additionally provides
polarity information.
"""
function srf_footprint_mask(
    func_pval::Function, raw::AbstractArray; alpha=0.01, fdr_c=:auto
)
    pval = func_pval(raw)
    C = if fdr_c == :auto
        benjamini_hochberg_constant(length(pval))
    else
        fdr_c
    end
    qval, qk = benjamini_hochberg_qvalue(pval; C)
    return srf_footprint_mask_from_qvalue(qval, invperm(qk[:]); alpha) #NOTE: invperm only takes a vector.
end

function srf_footprint_mask(raw::AbstractArray; kwargs...)
    return srf_footprint_mask(zscore_pvalue, raw; kwargs...)
end

#NOTE: when N > 50, the difference is less than 0.01
#NOTE: when N > 2000, the `sum` starts to be significantly slower than the approximation
@doc raw"""
    benjamini_hochberg_constant(N::Integer)

Compute the Benjamini–Hochberg correction constant for `N` tests.

Uses the exact harmonic sum when `N ≤ 50`, otherwise the approximation
`ln(N) + γ + 1/N`.
"""
function benjamini_hochberg_constant(N::Integer)
    return N > 50 ? log(N) + Base.MathConstants.eulergamma + 1 / N : sum(x -> 1 / x, 1:N)
end

@doc raw"""
    benjamini_hochberg_qvalue(pvalue; C=1) -> (q::Array, k::Array)

False discovery rate test for any list of pvalues.
This function will return the `q` value and the corresponding `k` index.

```math
P_k \leq \frac{k}{N} \alpha
```

```math
q_k = \frac{N}{k} P_k \leq \alpha
```

`k` is generated from `invperm(sortperm(...))`.
"""
function benjamini_hochberg_qvalue(pvalue::AbstractArray; C=1)
    N = length(pvalue)
    k = reshape(invperm(sortperm(pvalue[:])), size(pvalue)...)
    qval = pvalue .* (N * C) ./ k

    return (qval, k)
end

@doc raw"""
    zscore_pvalue(zscore::AbstractArray; two_tailed=true)

Convert z-scores to p-values under a standard normal distribution.

## Arguments
- `zscore`: array of z-scores.
- `two_tailed`: use two-tailed test (default `true`).
"""
function zscore_pvalue(zscore::AbstractArray; two_tailed=true)
    p = map(Base.Fix1(cdf, Normal()), -1 .* abs.(zscore))
    return two_tailed ? 2 .* p : p
end

@doc raw"""
    srf_footprint_mask_from_qvalue(qval, k; alpha=0.01)

Create a binary footprint mask from pre-computed q-values and their rank permutation
vector using the Benjamini–Hochberg procedure.

## Arguments
- `qval`: array of q-values.
- `k`: rank permutation vector (from `invperm(sortperm(...))`).
- `alpha`: FDR threshold (default `0.01`).
"""
function srf_footprint_mask_from_qvalue(qval, k; alpha=0.01)
    _cutoff = findlast(qval[k] .<= alpha)

    _mask = zeros(Bool, size(qval))
    if !isnothing(_cutoff)
        _mask[k[1:_cutoff]] .= true
    end

    return _mask
end
