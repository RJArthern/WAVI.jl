"""
    SolverParams(;
        n_iter_viscosity = 2,
        tol_picard = 1e-5,
        maxiter_picard = 30,
        tol_coarse = 1e-5,
        maxiter_coarse = 1000,
        levels = 3,
        wavelet_threshold = 10.0,
        nsmooth = 5,
        smoother_omega = 1.0,
        stencil_margin = 3)

Construct a WAVI.jl solver parameters object.
"""
@with_kw struct SolverParams{T <: Real, N <: Integer}
    n_iter_viscosity::N = 2;  @assert n_iter_viscosity ==2
    tol_picard::T = 1e-5
    maxiter_picard::N = 30
    tol_coarse::T = 1e-5
    maxiter_coarse::N = 1000
    levels::N = 3
    wavelet_threshold::T = 10.0
    nsmooth::N = 5
    smoother_omega::T = 1.0
    stencil_margin::N = 3
end