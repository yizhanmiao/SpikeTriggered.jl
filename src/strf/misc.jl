# misc functions
export floodfill!

@doc raw"""
	floodfill!(mat, initnode, target, replace; connectivity=4)

A [flood fill algorithm](https://rosettacode.org/wiki/Bitmap/Flood_fill) to
identify connected regions (values that equal to `target`) within a 2D matrix.

Valid elements will be replace with the `replace` value.

`connectivity` can be either `4` or `8`.

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

    queue = [initnode]

    seen = zeros(Bool, size(mat))
    seen[initnode] = true
    mat[initnode] = replace

    while !isempty(queue)
        current = pop!(queue)

        for δ in steps
            neighbor = CartesianIndex(current.I .+ δ)
            if inside_bounds(neighbor) && mat[neighbor] == target && !seen[neighbor]
                seen[neighbor] = true
                mat[neighbor] = replace
                push!(queue, neighbor)
            end
        end
    end

    return mat
end
