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
    ensure_colour_indices!(s, gu, gv)
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
        colour_indices=s.gs_colour_indices,
        resid_tmp=s.gs_resid,
        b_coarse=b_coarse,
        prolonged=s.prolonged,
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
            colour_indices,resid_tmp,b_coarse,prolonged = p

    # Multigrid smooth
    gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                           colour_indices=colour_indices, sweep_order=GS_SWEEP_ORDER,
                           smoother_omega=smoother_omega, resid=resid_tmp)

    # Multigrid restriction. resid_tmp is already b - A x after the last colour update.
    mul!(b_coarse, restrict, resid_tmp)

    abstol = tol_coarse * norm(b_coarse)

    # Multigrid solve for correction
    cg!(correction_coarse, op_coarse, b_coarse; abstol = abstol, maxiter = maxiter_coarse)

    # Multigrid prolongation
    mul!(prolonged, prolong, correction_coarse)
    @. x = x + prolonged

    # Multigrid smooth
    gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                           colour_indices=colour_indices, sweep_order=GS_SWEEP_ORDER_REV,
                           smoother_omega=smoother_omega, resid=resid_tmp)

    return x
end


"""
    get_op_diag(model, op)

The Gauss-Seidel weights are the force at each inner point if that point's
velocity is 1 and every other velocity is 0. They are written into stencil
scratch and overwritten on the next call.
"""
function get_op_diag(model::AbstractModel, op::LinearMap)
    @unpack gh, gu, gv, gc = model.fields
    s = stencil_scratch!(model)
    mi, ni = size(op)
    @assert mi == ni == gu.ni + gv.ni

    op_diag = s.op_diag
    T = eltype(op)
    dx_inv = one(T) / model.grid.dx
    dy_inv = one(T) / model.grid.dy
    D_h = reshape(gh.dneghηav[].diag, gh.nxh, gh.nyh)
    D_c = reshape(gc.dneghηav[].diag, gc.nxc, gc.nyc)
    D_u = reshape(gu.dnegβeff[].diag, gu.nxu, gu.nyu)
    D_v = reshape(gv.dnegβeff[].diag, gv.nxv, gv.nyv)
    D_imp = reshape(gh.dimplicit[].diag, gh.nxh, gh.nyh)

    if gu.ni > 0
        launch!(_op_diag_u!, view(op_diag, 1:gu.ni), s.gu_inner_indices,
                D_h, D_c, D_u, D_imp, gu.h, gu.mask, gc.mask, dx_inv, dy_inv;
                ndrange = gu.ni)
    end
    if gv.ni > 0
        launch!(_op_diag_v!, view(op_diag, (gu.ni + 1):ni), s.gv_inner_indices,
                D_h, D_c, D_v, D_imp, gv.h, gv.mask, gc.mask, dx_inv, dy_inv;
                ndrange = gv.ni)
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
    gauss_seidel_smoother!(x, op, b; iters, op_diag, colour_indices, sweep_order, smoother_omega, resid)

Apply smoother used in multigrid preconditioner.
`colour_indices[c]` is the packed list of degrees of freedom for colour `c`.
`resid` is workspace for `b - op * x`.
"""
function gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag,
                                colour_indices,
                                sweep_order,
                                smoother_omega=1.0,
                                resid)
    get_resid!(resid, x, op, b)
    for _ in 1:iters
        for j in sweep_order
            idx = colour_indices[j]
            isempty(idx) && continue
            @inbounds for k in idx
                x[k] += smoother_omega * resid[k] / op_diag[k]
            end
            get_resid!(resid, x, op, b)
        end
    end
    return x
end

function ensure_colour_indices!(s, gu, gv)
    if !s.gs_colours_filled
        s.gs_colour_indices = gs_colour_index_lists(gu, gv)
        s.gs_colours_filled = true
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

function ensure_multigrid_ops!(model::AbstractModel{T,N}, op::LinearMap{T}, s) where {T,N}
    @unpack wu, wv = model.fields
    n_wu = wu.n[]
    n_wv = wv.n[]
    n_coarse = n_wu + n_wv
    ni = size(op, 1)
    mg = s.mg_ops
    if mg.n_wu != n_wu || mg.n_wv != n_wv
        s.mg_ops = MultigridScratch{T}(;
            n_wu,
            n_wv,
            b_coarse = zeros(T, n_coarse),
            correction_coarse = zeros(T, n_coarse),
        )
        mg = s.mg_ops
    end
    restrict_fun! = get_restrict_fun(model)
    prolong_fun! = get_prolong_fun(model)
    restrict = LinearMap{T}(restrict_fun!, n_coarse, ni; issymmetric=false, ismutating=true, ishermitian=false, isposdef=false)
    prolong = LinearMap{T}(prolong_fun!, ni, n_coarse; issymmetric=false, ismutating=true, ishermitian=false, isposdef=false)
    op_coarse_fun! = get_op_coarse_fun(op, restrict, prolong, s.op_coarse_tmp1, s.op_coarse_tmp2)
    op_coarse = LinearMap{T}(op_coarse_fun!, n_coarse, n_coarse; issymmetric=true, ismutating=true, ishermitian=true, isposdef=true)
    return restrict, prolong, op_coarse, mg.b_coarse, mg.correction_coarse
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

function get_op_coarse_fun(
    op::LinearMap{T},
    restrict::LinearMap{T},
    prolong::LinearMap{T},
    tmp1::AbstractVector{T},
    tmp2::AbstractVector{T},
) where {T}
    mi, ni = size(op)
    @assert mi == ni == length(tmp1) == length(tmp2)

    function op_coarse_fun!(out, in)
        mul!(tmp1, prolong, in)
        mul!(tmp2, op, tmp1)
        mul!(out, restrict, tmp2)
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

