# misc functions

@doc raw"""
    floodfill!(mat, initnode, target, replace; connectivity=4)

Flood-fill connected region starting from `initnode`, replacing `target` values
with `replace`.

## Arguments
- `mat::AbstractMatrix`: matrix to modify in-place.
- `initnode::CartesianIndex{2}`: starting position.
- `target`: value to match for flood filling.
- `replace`: replacement value.
- `connectivity`: `4` or `8` (default `4`).

Returns the modified `mat`.
"""
function floodfill!(
    mat::AbstractMatrix, initnode::CartesianIndex{2}, target, replace; connectivity=4
)
    mat[initnode] == target || return mat

    (w, h) = size(mat)
    inside_bounds(i) = (1 <= i.I[1] <= w) && (1 <= i.I[2] <= h)

    steps = if connectivity == 4
        [(-1, 0), (1, 0), (0, -1), (0, 1)]   # 4-connectivity
    elseif connectivity == 8
        [(i, j) for i in -1:1, j in -1:1 if !(i==0 && j==0)]  # 8-connectivity
    end

    stack = [initnode]

    seen = zeros(Bool, size(mat))
    seen[initnode] = true
    mat[initnode] = replace

    while !isempty(stack)
        current = pop!(stack)

        for δ in steps
            neighbor = CartesianIndex(current.I .+ δ)
            if inside_bounds(neighbor) && mat[neighbor] == target && !seen[neighbor]
                seen[neighbor] = true
                mat[neighbor] = replace
                push!(stack, neighbor)
            end
        end
    end

    return mat
end
