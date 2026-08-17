export apply_preconditioning!, gauss_seidel_smoother!, get_correction_coarse,
    get_multigrid_ops, get_op_coarse_fun, get_op_diag, get_preconditioner,
    set_correction_coarse!

using IterativeSolvers

using WAVI.Utilities

# Red-black colours: u on 1–2, v on 3–4. Post-smooth uses the reverse order.
const GS_SWEEP_ORDER = [1, 3, 2, 4]
const GS_SWEEP_ORDER_REV = [4, 2, 3, 1]

"""
    get_preconditioner(model::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

 Get precondtioner. Details are stored in a struct of type Preconditioner that can be passed to a solver.

"""
function get_preconditioner(model::AbstractModel{T,N}, op::LinearMap{T}) where {T,N}
    @unpack gu,gv=model.fields
    @unpack solver_params=model
    s = stencil_scratch!(model)

    mi,ni = size(op)
    @assert mi == ni == gu.ni + gv.ni

    restrict, prolong, op_coarse, b_coarse, correction_coarse = ensure_multigrid_ops!(model, op, s)
    #starting guess for the multigrid coarse correction is cached by the model
    fill_correction_coarse!(model, correction_coarse)
    op_diag = get_op_diag(model, op)

    #Four colour Jacobi preconditioner. Red-Black checkerboard Jacobi for each velocity component.
    O=typeof(op)
    C=typeof(op_coarse)
    R=typeof(restrict)
    P=typeof(prolong)
    p=Preconditioner{T,N,O,C,R,P}(
        op=op,
        restrict=restrict,
        prolong=prolong,
        op_coarse=op_coarse,
        op_diag=op_diag,
        nsmooth=solver_params.nsmooth,
        tol_coarse=solver_params.tol_coarse,
        correction_coarse=correction_coarse,
        maxiter_coarse=solver_params.maxiter_coarse,
        smoother_omega=solver_params.smoother_omega,
        colour_indices=s.gs_colour_indices[],
        resid_tmp=s.gs_resid,
        b_coarse=b_coarse,
        prolonged=s.prolonged,
        gs_increment=s.gs_increment,
        gs_applied_increment=s.gs_applied_increment,
        apply_colour=s.apply_colour_op!,
    )

    return p
end

function get_correction_coarse(p::AbstractPreconditioner)
    return p.correction_coarse
end

"""
    apply_preconditioning!(x, p, b)

Apply wavelet-based multigrid preconditioner using information stored in p.
"""
function apply_preconditioning!(x, p, b)
    @unpack op,op_diag,nsmooth,smoother_omega,restrict,
            prolong,op_coarse,correction_coarse,tol_coarse,maxiter_coarse,
            colour_indices,resid_tmp,b_coarse,prolonged,
            gs_increment,gs_applied_increment,apply_colour = p

    # Preconditioner apply starts from x = 0, so the residual is b.
    fill!(x, zero(eltype(x)))
    copyto!(resid_tmp, b)
    gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                           colour_indices=colour_indices, sweep_order=GS_SWEEP_ORDER,
                           smoother_omega=smoother_omega, resid=resid_tmp,
                           increment=gs_increment, applied_increment=gs_applied_increment, apply_colour=apply_colour,
                           resid_is_current=true)

    # Multigrid restriction. resid_tmp is already b - A x after the last colour update.
    mul!(b_coarse, restrict, resid_tmp)

    abstol = tol_coarse * norm(b_coarse)

    # Multigrid solve for correction
    cg!(correction_coarse, op_coarse, b_coarse; abstol = abstol, maxiter = maxiter_coarse)

    # Multigrid prolongation
    mul!(prolonged, prolong, correction_coarse)
    @. x = x + prolonged

    # x changed by the coarse correction; rebuild the opening residual.
    gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                           colour_indices=colour_indices, sweep_order=GS_SWEEP_ORDER_REV,
                           smoother_omega=smoother_omega, resid=resid_tmp,
                           increment=gs_increment, applied_increment=gs_applied_increment, apply_colour=apply_colour,
                           resid_is_current=false)

    return x
end


"""
    get_op_diag(wavi::AbstractModel,op::LinearMap)

 Get diagonal of operator for use in preconditioner.
 The returned vector is stencil scratch and is overwritten on the next call.

"""
function get_op_diag(model::AbstractModel,op::LinearMap)
    @unpack gu,gv=model.fields
    s = stencil_scratch!(model)
    mi,ni = size(op)
    @assert mi == ni == gu.ni + gv.ni

    op_diag = s.op_diag
    probe = s.diag_probe
    tmp = s.diag_tmp
    ensure_colour_indices!(s, gu, gv, model.solver_params.stencil_margin)
    fill!(op_diag, zero(eltype(op)))
    for idx in s.diag_colour_indices[]
        isempty(idx) && continue
        fill!(probe, zero(eltype(probe)))
        @inbounds for k in idx
            probe[k] = one(eltype(probe))
        end
        mul!(tmp, op, probe)
        @inbounds for k in idx
            op_diag[k] = tmp[k]
        end
    end
    return op_diag
end


"""
    LinearAlgebra.ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}

Overload LinearAlgebra.ldiv! function so that the bespoke preconditioner is deployed in calls to the
conjugate gradient method if p has type  <: AbstractPreconditioner.

"""
function ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T}, b::AbstractVecOrMat{T}) where {T}
    apply_preconditioning!(x, p, b)
end


"""
    gauss_seidel_smoother!(x, op, b; iters, op_diag, colour_indices, sweep_order,
                           smoother_omega, resid, increment, applied_increment,
                           apply_colour, resid_is_current)

Apply smoother used in multigrid preconditioner.
`colour_indices[c]` is the packed list of degrees of freedom for colour `c`.
`resid` is workspace for `b - op * x`. After each colour, only that colour's
increment is applied through `apply_colour`, then
`resid = resid - A * increment`.
If `resid_is_current` is true, `resid` is already `b - op * x` and the opening
matvec is skipped.
"""
function gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag,
                                colour_indices,
                                sweep_order,
                                smoother_omega=1.0,
                                resid,
                                increment,
                                applied_increment,
                                apply_colour,
                                resid_is_current=false)
    if !resid_is_current
        get_resid!(resid, x, op, b)
    end
    for _ in 1:iters
        for j in sweep_order
            idx = colour_indices[j]
            isempty(idx) && continue
            launch!(_gs_colour_saxpy!, x, increment, resid, op_diag, idx, smoother_omega;
                    ndrange = length(idx))
            apply_colour(applied_increment, increment, idx)
            launch!(_gs_resid_sub!, resid, applied_increment; ndrange = length(resid))
        end
    end
    return x
end

function ensure_colour_indices!(s, gu, gv, sm)
    if s.gs_colour_indices[] === nothing
        s.gs_colour_indices[] = gs_colour_index_lists(gu, gv)
    end
    if s.diag_colour_indices[] === nothing
        s.diag_colour_indices[] = diag_colour_index_lists(gu, gv, sm)
    end
    return nothing
end

"""
Packed red-black colours for Gauss-Seidel: u on 1–2, v on 3–4.
Order is column-major over `mask_inner`, matching packed `samp_inner`.
"""
function gs_colour_index_lists(gu, gv)
    lists = [Int[] for _ in 1:4]
    k = 0
    for j in 1:gu.nyu, i in 1:gu.nxu
        if gu.mask_inner[i, j]
            k += 1
            push!(lists[1 + mod(i - j, 2)], k)
        end
    end
    for j in 1:gv.nyv, i in 1:gv.nxv
        if gv.mask_inner[i, j]
            k += 1
            push!(lists[3 + mod(i - j, 2)], k)
        end
    end
    return lists
end

"""
Packed colours for Jacobi diagonal probes, spaced by `stencil_margin`
so neighbouring degrees of freedom are not in the same probe.
"""
function diag_colour_index_lists(gu, gv, sm)
    lists = [Int[] for _ in 1:(2 * sm^2)]
    k = 0
    for j in 1:gu.nyu, i in 1:gu.nxu
        if gu.mask_inner[i, j]
            k += 1
            push!(lists[1 + mod(i - 1, sm) + sm * mod(j - 1, sm)], k)
        end
    end
    for j in 1:gv.nyv, i in 1:gv.nxv
        if gv.mask_inner[i, j]
            k += 1
            push!(lists[1 + sm^2 + mod(i - 1, sm) + sm * mod(j - 1, sm)], k)
        end
    end
    return lists
end

function ensure_multigrid_ops!(model::AbstractModel{T,N}, op::LinearMap{T}, s) where {T,N}
    @unpack wu, wv = model.fields
    n_wu = wu.n[]
    n_wv = wv.n[]
    n_coarse = n_wu + n_wv
    mg = s.mg_ops[]
    if mg !== nothing && mg.n_wu == n_wu && mg.n_wv == n_wv
        return mg.restrict, mg.prolong, mg.op_coarse, mg.b_coarse, mg.correction_coarse
    end

    ni = size(op, 1)
    restrict_fun! = get_restrict_fun(model)
    prolong_fun! = get_prolong_fun(model)
    restrict = LinearMap{T}(restrict_fun!, n_coarse, ni; issymmetric=false, ismutating=true, ishermitian=false, isposdef=false)
    prolong = LinearMap{T}(prolong_fun!, ni, n_coarse; issymmetric=false, ismutating=true, ishermitian=false, isposdef=false)
    op_coarse_fun! = get_op_coarse_fun(op, restrict, prolong)
    op_coarse = LinearMap{T}(op_coarse_fun!, n_coarse, n_coarse; issymmetric=true, ismutating=true, ishermitian=true, isposdef=true)

    s.mg_ops[] = (
        n_wu = n_wu,
        n_wv = n_wv,
        restrict = restrict,
        prolong = prolong,
        op_coarse = op_coarse,
        b_coarse = zeros(T, n_coarse),
        correction_coarse = zeros(T, n_coarse),
    )
    mg = s.mg_ops[]
    return mg.restrict, mg.prolong, mg.op_coarse, mg.b_coarse, mg.correction_coarse
end

function fill_correction_coarse!(model::AbstractModel, correction_coarse)
    @unpack wu, wv = model.fields
    n_coarse = wu.n[] + wv.n[]
    correction_coarse[1:wu.n[]] .= wu.correction_coarse[]
    correction_coarse[(wu.n[] + 1):n_coarse] .= wv.correction_coarse[]
    return correction_coarse
end

function get_multigrid_ops(model::AbstractModel{T,N}, op::LinearMap{T}) where {T,N}
    s = stencil_scratch!(model)
    restrict, prolong, op_coarse, _, _ = ensure_multigrid_ops!(model, op, s)
    return restrict, prolong, op_coarse
end

function get_op_coarse_fun(op::LinearMap{T},restrict::LinearMap{T},prolong::LinearMap{T}) where {T}

     mi,ni = size(op)
     @assert mi == ni 

     tmp1 :: Vector{T} = zeros(ni)
     tmp2 :: Vector{T} = zeros(ni)
     function op_coarse_fun!(out,in)
@!        tmp1 = prolong * in
@!        tmp2 = op * tmp1
@!        out = restrict * tmp2
          return out
     end
     return op_coarse_fun!
end

function set_correction_coarse!(model::AbstractModel{T,N},correction_coarse::AbstractVector{T}) where {T,N}
    @unpack wu,wv = model.fields

    n_coarse = wu.n[] + wv.n[]

    wu.correction_coarse[] .= correction_coarse[ 1 : wu.n[] ]  
    wv.correction_coarse[] .= correction_coarse[ (wu.n[]+1) : n_coarse]

    return nothing
end

