
struct SparseNoise{T}
    ids::Vector{Int16}
    row::Vector{Int16}
    col::Vector{Int16}
    sign::Vector{T}
    gridsize::Integer

    function SparseNoise(snra::Vector{T}, gridsize::Integer) where {T<:Integer}
        _snra = parse_snra(snra, gridsize)
        new{eltype(_snra.sign)}(_snra.ids, _snra.row, _snra.col, _snra.sign, gridsize)
    end

    function SparseNoise(snf_file::String, args...; kwargs...)
        _snra = load_snra(snf_file)
        SparseNoise(_snra, args...; kwargs...)
    end
end

@doc raw"""
    load_snra(snf_path::String) -> Vector{Int16}

read `snf` sparse noise pseudorandom files

## Synopsis

```c
struct SNF {
    int16_t n_order;
    int32_t base_poly; // not sure what it represents
    int16_t* snra; // the length would be 2^n_order
};
```
"""
function load_snra(snf_path::String)
    @assert isfile(snf_path) && endswith(snf_path, ".snf") "input must be valid SNF file path."

    fid = open(snf_path)
    n_order = read(fid, Int16)
    base_poly = read(fid, Int32) #XXX: not sure what bose_ploy represents.
    snra = collect(reinterpret(Int16, read(fid)))
    close(fid)

    length(snra) == 2 ^ n_order ||
        @warn "snf file data length not match with the record" length(snra) 2^n_order

    snra
end

@doc raw"""
    parse_snra(snra::Array{T}, grid_size::Integer) -> NamedTuple(:ids, :row, :col, :sign)

parse the snra array to sparse noise image info.
"""
function parse_snra(snra::Array{T}, grid_size::Integer) where {T<:Integer}
    _aftermod = mod.(snra .- 1, grid_size * grid_size * 2)
    _row = _aftermod .÷ (grid_size * 2) .+ 1
    _col_n = mod.(_aftermod, grid_size * 2)
    _col = _col_n .÷ 2 .+ 1
    _sign = 1 .- mod.(_col_n, 2) .* 2

    return (ids=T.(_aftermod), row=T.(_row), col=T.(_col), sign=T.(_sign))
end

function Base.collect(sn::SparseNoise{T}) where {T}
    N = length(sn.ids)
    output = zeros(T, N, sn.gridsize * sn.gridsize)
    for idx in 1:N
        # @inbounds output[idx, (sn.row[idx]-1)*sn.gridsize+sn.col[idx]] = sn.sign[idx]
        @inbounds output[idx, (sn.col[idx] - 1) * sn.gridsize + sn.row[idx]] = sn.sign[idx]
    end
    return output
end
