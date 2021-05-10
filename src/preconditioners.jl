"""
    get_preconditioner(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

 Get precondtioner. Details are stored in a struct of type Preconditioner that can be passed to a solver.

"""
function get_preconditioner(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

    @unpack gu,gv,wu,wv,params,solver_params=model

    m,n=size(op)
    @assert m == n == gu.n + gv.n

    n = gu.n + gv.n
    n_coarse = wu.n[] + wv.n[]

    restrict_fun(x) = restrictvec(model,x)
    prolong_fun(x) = prolongvec(model,x)

    restrict=LinearMap{T}(restrict_fun,n_coarse,n;issymmetric=false,ismutating=false,ishermitian=false,isposdef=false)
    prolong=LinearMap{T}(prolong_fun,n,n_coarse;issymmetric=false,ismutating=false,ishermitian=false,isposdef=false)

    op_diag=get_op_diag(model,op)

    #Four colour Jacobi preconditioner. Red-Black checkerboard Jacobi for each velocity component.
    sweep=[[1 .+ mod(i-j,2) for i=1:gu.nx, j=1:gu.ny][gu.mask];[3 .+ mod(i-j,2) for i=1:gv.nx, j=1:gv.ny][gv.mask]]
    sweep_order=[1,3,2,4]

    p=Preconditioner{T,N}(op=op, restrict=restrict, prolong=prolong,sweep=sweep, sweep_order=sweep_order,
            op_diag=op_diag, nsmooth=solver_params.nsmooth, tol_coarse = solver_params.tol_coarse,
            maxiter_coarse = solver_params.maxiter_coarse, smoother_omega=solver_params.smoother_omega)

    return p
end

"""
    precondition!(x, p, b)

Apply wavelet-based multigrid preconditioner using information stored in p.
"""
function precondition!(x, p, b)
    @unpack op,op_diag,nsmooth,sweep,sweep_order,smoother_omega,restrict,
            prolong,op_coarse,correction_coarse,tol_coarse,maxiter_coarse = p


    n=size(op,1)

    # Multigrid smooth
    x .= gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                                sweep=sweep, sweep_order=sweep_order, smoother_omega = smoother_omega)

    resid=b-op*x;

    # Multigrid restriction
    b_coarse=restrict*resid

    # Multigrid solve for correction
    cg!(correction_coarse, op_coarse, b_coarse; reltol = tol_coarse, maxiter = maxiter_coarse)

    # Multigrid prolongation
    x .= x .+ prolong*correction_coarse

    # Multigrid smooth
    x .= gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                                sweep=sweep, sweep_order=reverse(sweep_order), smoother_omega = smoother_omega)

    return x
end


"""
    get_op_diag(wavi::AbstractModel,op::LinearMap)

 Get diagonal of operator for use in preconditioner.

"""
function get_op_diag(wavi::AbstractModel,op::LinearMap)
    @unpack gu,gv,params,solver_params=wavi
    m,n=size(op)
    @assert m == n == gu.n + gv.n
    op_diag=zeros(eltype(op),n)
    sm=solver_params.stencil_margin
    sweep=[[1+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gu.nx, j=1:gu.ny][gu.mask];
           [1+sm^2+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gv.nx, j=1:gv.ny][gv.mask] ]
    e=zeros(Bool,n)
    for i = unique(sweep)
        e[sweep .== i] .= true
        op_diag[e] .= (op*e)[e]
        e[sweep .== i] .= false
    end
    return op_diag
end


"""
    LinearAlgebra.ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}

Overload LinearAlgebra.ldiv! function so that the bespoke preconditioner is deployed in calls to the
conjugate gradient method if p has type  <: AbstractPreconditioner.

"""
function ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}
    precondition!(x, p, b)
end


"""
    gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag=diag(sparse(op)),
                                sweep=(1:size(op,1)),
                                sweep_order=unique(sweep),
                                smoother_omega=1.0)
Apply smoother used in multigrid preconditioner.
"""
function gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag=diag(sparse(op)),
                                sweep=(1:size(op,1)),
                                sweep_order=unique(sweep),
                                smoother_omega=1.0)
    resid=b-op*x
    for i = 1:iters
        for j = sweep_order
              idx = sweep .== j
              x[idx] .= x[idx] .+ smoother_omega .* resid[idx]./op_diag[idx]
              resid .= b .- op*x
        end
    end
    return x
end