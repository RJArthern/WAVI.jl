"""
    get_preconditioner(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

 Get precondtioner. Details are stored in a struct of type Preconditioner that can be passed to a solver.

"""
function get_preconditioner(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

    @unpack gu,gv,wu,wv=model.fields
    @unpack params,solver_params=model

    m,n=size(op)
    @assert m == n == gu.n + gv.n

    n = gu.n + gv.n
    n_coarse = wu.n[] + wv.n[]

    #starting guess for the multigrid coarse correction is cached by the model
    correction_coarse=get_correction_coarse(model)

    restrict,prolong,op_coarse = get_multigrid_ops(model,op)

    op_diag=get_op_diag(model,op)

    #Four colour Jacobi preconditioner. Red-Black checkerboard Jacobi for each velocity component.
    sweep=[[1 .+ mod(i-j,2) for i=1:gu.nxu, j=1:gu.nyu][gu.mask];[3 .+ mod(i-j,2) for i=1:gv.nxv, j=1:gv.nyv][gv.mask]]
    sweep_order=[1,3,2,4]

    O=typeof(op)
    C=typeof(op_coarse)
    R=typeof(restrict)
    P=typeof(prolong)
    p=Preconditioner{T,N,O,C,R,P}(op=op, restrict=restrict, prolong=prolong, op_coarse = op_coarse, sweep=sweep, sweep_order=sweep_order,
            op_diag=op_diag, nsmooth=solver_params.nsmooth, tol_coarse = solver_params.tol_coarse, correction_coarse = correction_coarse,
            maxiter_coarse = solver_params.maxiter_coarse, smoother_omega=solver_params.smoother_omega)

    return p
end

function get_correction_coarse(p::AbstractPreconditioner)
    return p.correction_coarse
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

    resid=get_resid(x,op,b)

    # Multigrid restriction
    b_coarse=restrict*resid

    abstol = tol_coarse*norm(b_coarse)

    # Multigrid solve for correction
    cg!(correction_coarse, op_coarse, b_coarse; abstol = abstol, maxiter = maxiter_coarse)

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
function get_op_diag(model::AbstractModel,op::LinearMap)
    @unpack gu,gv=model.fields
    @unpack params,solver_params=model
    m,n=size(op)
    @assert m == n == gu.n + gv.n
    op_diag=zeros(eltype(op),n)
    ope_tmp=zeros(eltype(op),n)
    sm=solver_params.stencil_margin
    sweep=[[1+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gu.nxu, j=1:gu.nyu][gu.mask];
           [1+sm^2+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gv.nxv, j=1:gv.nyv][gv.mask] ]
    e=zeros(Bool,n)
    for i = unique(sweep)
        e .= sweep .== i
        mul!(ope_tmp,op,e)
        op_diag[e] .= ope_tmp[e]
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
    resid=get_resid(x,op,b)
    
    idx = zeros(Bool,size(sweep))
    for i = 1:iters
        for j = sweep_order
              idx .= sweep .== j
              x[idx] .= x[idx] .+ smoother_omega .* resid[idx]./op_diag[idx]
              get_resid!(resid,x,op,b)
        end
    end
    return x
end


function get_multigrid_ops(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}
    @unpack gu,gv,wu,wv=model.fields

    m,n=size(op)
    @assert m == n == gu.n + gv.n

    n = gu.n + gv.n
    n_coarse = wu.n[] + wv.n[]

    restrict_fun! = get_restrict_fun(model)
    prolong_fun! = get_prolong_fun(model)

    restrict=LinearMap{T}(restrict_fun!,n_coarse,n;issymmetric=false,ismutating=true,ishermitian=false,isposdef=false)
    prolong=LinearMap{T}(prolong_fun!,n,n_coarse;issymmetric=false,ismutating=true,ishermitian=false,isposdef=false)

    op_coarse_fun! = get_op_coarse_fun(op,restrict,prolong)

    op_coarse=LinearMap{T}(op_coarse_fun!,n_coarse,n_coarse;issymmetric=true,ismutating=true,ishermitian=true,isposdef=true)

    return restrict,prolong,op_coarse

end

function get_op_coarse_fun(op::LinearMap{T},restrict::LinearMap{T},prolong::LinearMap{T}) where {T}

     m,n = size(op)
     @assert m == n 

     n_coarse = size(restrict,1)

     tmp1 :: Vector{T} = zeros(n)
     tmp2 :: Vector{T} = zeros(n)
     out :: Vector{T} = zeros(n_coarse)
     function op_coarse_fun!(out,in)
@!        tmp1 = prolong * in
@!        tmp2 = op * tmp1
@!        out = restrict * tmp2
          return out
     end
     return op_coarse_fun!
end

function get_correction_coarse(model::AbstractModel{T,N}) where {T,N}
    @unpack wu,wv = model.fields

    n_coarse = wu.n[] + wv.n[]

    correction_coarse=zeros(T,n_coarse)
    correction_coarse[1:wu.n[]] .= wu.correction_coarse[]
    correction_coarse[wu.n[]+1:n_coarse] .= wv.correction_coarse[]

    return correction_coarse
end

function set_correction_coarse!(model::AbstractModel{T,N},correction_coarse::AbstractVector{T}) where {T,N}
    @unpack wu,wv = model.fields

    n_coarse = wu.n[] + wv.n[]

    wu.correction_coarse[] .= correction_coarse[ 1 : wu.n[] ]  
    wv.correction_coarse[] .= correction_coarse[ (wu.n[]+1) : n_coarse]

    return nothing
end

