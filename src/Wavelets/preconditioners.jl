#Struct to hold information about wavelet-based multigrid preconditioner.
@with_kw struct Preconditioner{T <: Real, N <: Integer} <: AbstractPreconditioner{T,N}
    op::LinearMap{T}
    op_diag::Vector{T} = diag(sparse(op))
    nsmooth::N = 5
    sweep::Vector{N}
    sweep_order::Vector{N} = unique(sweep)
    smoother_omega::T = 1.0
    restrict::LinearMap{T}
    prolong::LinearMap{T}
    op_coarse::LinearMap{T} = restrict*op*prolong
    correction_coarse::Vector{T} = zeros(T,size(op_coarse,2))
    tol_coarse::T = 1e-7
    maxiter_coarse::N = 1000
end
