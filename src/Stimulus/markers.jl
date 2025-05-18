#TODO: docs
function process_first_of_paired_markers(m::Vector{T}) where {T}
    @assert length(m) % 2 == 0 "markers have to be paired!"
    m[1:2:(end - 1)]
end

function remove_invalid_serial_timestamps!(
    markers::Vector{T}; t_interval=1/35, atol=1e-1, interior=false
) where {T}
    Δt = diff(markers)
    _bad_marker = Int64[]

    # check the beginning
    _found_start = false
    for idx in 1:(length(markers) - 1)
        if !isapprox(Δt[idx], t_interval; atol)
            push!(_bad_marker, idx)
            _found_start = true
        elseif _found_start && !interior
            break
        end
    end

    if _found_start && !interior
        _found_end = false
        for idx in (length(markers) - 1):(-1):(maximum(_bad_marker) + 1)
            if !isapprox(Δt[idx], t_interval; atol)
                push!(_bad_marker, idx+1)
                _found_end = true
            elseif _found_end
                break
            end
        end
    end

    _serial_idx = cumsum(ones(length(_bad_marker)))
    for idx in sort(_bad_marker) .- _serial_idx .+ 1
        popat!(markers, Int64(idx))
    end

    markers, _bad_marker
end

function remove_invalid_paired_timestamps!(
    markers::Vector{T}; t_trial=2, t_interval=5, atol=1e-1, n_total
) where {T}
    Δt = diff(markers)
    _times = T[t_trial, t_interval]
    _status = 1; # status denote the position of _times
    _bad_marker = Int64[]

    for idx in 1:length(_times)
        if isapprox(Δt[idx], _times[_status]; atol)
            _status = _status % 2 + 1
        else
            push!(_bad_marker, idx)
        end
    end

    _serial_idx = cumsum(ones(length(_bad_marker)))
    for idx in sort(_bad_marker) .- _serial_idx .+ 1
        popat!(markers, Int64(idx))
    end

    if length(markers) % 2 == 0  #NOTE: paired
        markers
    else #NOTE: remove the last one to make it paired.
        pop!(markers)
        markers
    end

    if !isnothing(n_total) && length(markers) > 2*n_total
        markers = markers[1:(2 * n_total)]
    end

    markers, _bad_marker
end

function remove_invalid_timestamps(
    markers::Vector{T}; mtype::Symbol, debug=false, kwargs...
) where {T}
    tmp = copy(markers)
    output, popped_indice = if mtype == :serial
        remove_invalid_serial_timestamps!(tmp; kwargs...)
    elseif mtype == :paired
        remove_invalid_paired_timestamps!(tmp; kwargs...)
    else
        @error("unknown mtype: $(mtype)")
    end

    if debug
        output, popped_indice
    else
        output
    end
end

function get_trial_markers(raw_markers::AbstractVector; skip=1, step=2, length=20)
    raw_markers[range(skip+1; step, length)]
end
