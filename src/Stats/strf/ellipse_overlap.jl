# module EEOver

import Libdl

libEEOver = C_NULL
func_eeover = fill(C_NULL, 3)

function __init_eeover__()
    global libEEOver
    global func_eeover
    try
        libEEOver = Libdl.dlopen(joinpath(@__DIR__, "..", "..", "..", "deps", "eeover-c", "libeeover"))
        func_eeover[1] = Libdl.dlsym_e(libEEOver, :ellipse_ellipse_overlap_gsl)
        func_eeover[2] = Libdl.dlsym_e(libEEOver, :ellipse_ellipse_overlap_netlibs)
        func_eeover[3] = Libdl.dlsym_e(libEEOver, :ellipse_ellipse_overlap_gems)
    catch
        @warn("EEOver library not available!")
    end
end

# struct EEOverEllipse
#     axis_major::Cdouble
#     axis_minor::Cdouble
#     center_x::Cdouble
#     center_y::Cdouble
#     rotation::Cdouble
# end

const EEOverEllipse = GaussianEllipse;

struct EEOverResult
    overlapArea::Cdouble
    roots::Vector{NTuple{2, Cdouble}}
    returnCode::Cint
end

function eeover(e1::EEOverEllipse, e2::EEOverEllipse; solver=:gsl)
    root_x = Cdouble[0,0,0,0]
    root_y = Cdouble[0,0,0,0]
    NROOTS = Cint[0]
    rtnCode = Cint[0]

    func, choice = if solver == :gsl
        func_eeover[1], 1
    elseif solver == :optimal
        func_eeover[1], 2
    elseif solver == :tom
        func_eeover[2], nothing
    elseif solver == :gem
        func_eeover[3], nothing
    else
        @error("unknown solver option ($solver)")
    end

    rez = if choice isa Nothing
        @ccall $func(
            e1.rotation::Cdouble,
            e1.axis_major::Cdouble, e1.axis_minor::Cdouble,
            e1.center_x::Cdouble, e1.center_y::Cdouble,
            e2.rotation::Cdouble,
            e2.axis_major::Cdouble, e2.axis_minor::Cdouble,
            e2.center_x::Cdouble, e2.center_y::Cdouble,
            pointer(root_x)::Ptr{Cdouble}, pointer(root_y)::Ptr{Cdouble},
            pointer(NROOTS)::Ptr{Cint}, pointer(rtnCode)::Ptr{Cint},
        )::Cdouble
    else
        @ccall $func(
            e1.rotation::Cdouble,
            e1.axis_major::Cdouble, e1.axis_minor::Cdouble,
            e1.center_x::Cdouble, e1.center_y::Cdouble,
            e2.rotation::Cdouble,
            e2.axis_major::Cdouble, e2.axis_minor::Cdouble,
            e2.center_x::Cdouble, e2.center_y::Cdouble,
            pointer(root_x)::Ptr{Cdouble}, pointer(root_y)::Ptr{Cdouble},
            pointer(NROOTS)::Ptr{Cint}, pointer(rtnCode)::Ptr{Cint},
            choice::Cint
        )::Cdouble
    end

    n_roots = min(4, NROOTS[1])
    roots = [(root_x[idx], root_y[idx]) for idx in 1:n_roots]
    EEOverResult(rez, roots, rtnCode[1])
end

function ellipse_overlap_index(e1::GaussianEllipse, e2::GaussianEllipse; kwargs...)
    _overlap_area = eeover(e1, e2; kwargs...)
    2 * _overlap_area / (area(e1) + area(e2))
end

# end