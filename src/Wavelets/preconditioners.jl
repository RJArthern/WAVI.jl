#Struct to hold information about wavelet-based multigrid preconditioner.
@with_kw struct Preconditioner{T <: Real, N <: Integer, O <:MapOrMatrix{T}, C <: MapOrMatrix{T}, 
         R <: MapOrMatrix{T}, P <: MapOrMatrix{T}} <: AbstractPreconditioner{T,N}
    op::O
    op_diag::Vector{T} = diag(sparse(op))
    nsmooth::N = 5
    sweep::Vector{N}
    sweep_order::Vector{N} = unique(sweep)
    smoother_omega::T = 1.0
    restrict::R
    prolong::P
    op_coarse::C = restrict*op*prolong
    correction_coarse::Vector{T} = zeros(T,size(op_coarse,2))
    tol_coarse::T = 1e-7
    maxiter_coarse::N = 1000
end
