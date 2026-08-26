module Utilities

using LinearAlgebra
using Parameters
using SparseArrays

using WAVI: AbstractModel
using KernelAbstractions: KernelAbstractions as KA, @kernel, @index
using WAVI.Stencils

export get_op_fun, get_restrict_fun, get_prolong_fun, pos_fraction, pos_fraction!, mismip_plus_bed,
    get_glx, glen_b, fill_glen_b!, get_u_mask, get_v_mask, get_c_mask, clip, get_resid, get_resid!,
    icedraft, height_above_floatation, volume_above_floatation, spI, ∂1d, c, χ,
    stencil_scratch!, StencilScratch, apply_momentum_op!, copy_like, copy_onto!, zeros_like, _host

#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))

"""
    fill_index_map!(idx, mask) -> n

Build a 2D lookup from a true/false mask of the same size.
True cells are numbered 1, 2, 3, ... down each column, then left to right
(the position in the short packed vector). False cells get `0`.
Returns how many cells were true.
"""
function fill_index_map!(idx::AbstractMatrix{<:Integer}, mask::AbstractMatrix{Bool})
    if _is_host_array(idx) && _is_host_array(mask)
        return _fill_index_map_serial!(idx, mask)
    end
    idx_h = Array(idx)
    n = _fill_index_map_serial!(idx_h, Array(mask))
    copyto!(idx, idx_h)
    return n
end

_is_host_array(x) = x isa Array || x isa BitArray

function _fill_index_map_serial!(idx, mask)
    k = 0
    @inbounds for j in axes(mask, 2), i in axes(mask, 1)
        if mask[i, j]
            k += 1
            idx[i, j] = k
        else
            idx[i, j] = 0
        end
    end
    return k
end

"""
    inner_index_map(mask_inner) -> Matrix{Int}

The velocity solver stores one value per inner point rather than per grid
point. Inner points are ice that is free to move, so not ocean and not a
fixed boundary; `mask_inner` is true there.

This returns a 2D array of the same size as the grid. Where the point is
inner, `idx[i, j]` is its position `k` in the solver vector, and the
velocity is `u[k]`. Where it is not, `idx[i, j]` is `0`.

The list `inner_indices` runs the other way: `inner_indices[k]` is the 2D
location of the k-th inner point. Both number cells down each column, then
left to right, matching `findall(vec(mask_inner))`.
"""
function inner_index_map(mask_inner::AbstractMatrix{Bool})
    idx = zeros(Int, size(mask_inner))
    fill_index_map!(idx, mask_inner)
    return idx
end


"""
Coarse-grid work vectors for the wavelet multigrid cycle.
Resized when the number of kept wavelet coefficients changes.
"""
Base.@kwdef struct MultigridScratch{T <: Real}
    n_wu::Int
    n_wv::Int
    b_coarse::AbstractVector{T}
    correction_coarse::AbstractVector{T}
end

"""
Workspace for the momentum operator, Picard stencils, Haar RAP, and Gauss-Seidel.

Allocated once per `GridField` and reused. Colour lists are filled on first
preconditioner call. `op_coarse_tmp1` and `op_coarse_tmp2` are the RAP
coarse-operator work vectors. Coarse correction vectors sit in `mg_ops` and
are resized when the number of kept wavelets changes.
"""
Base.@kwdef mutable struct StencilScratch{T <: Real}
    gu_inner_indices::AbstractVector{Int}
    gv_inner_indices::AbstractVector{Int}
    surf_crop::AbstractArray{T, 2}
    ones_crop::AbstractArray{T, 2}
    tmpu::AbstractArray{T, 2}
    tmpv::AbstractArray{T, 2}
    tmpui::AbstractVector{T}
    tmpvi::AbstractVector{T}
    u_crop::AbstractArray{T, 2}
    v_crop::AbstractArray{T, 2}
    u_h::AbstractArray{T, 2}
    v_h::AbstractArray{T, 2}
    dudx::AbstractArray{T, 2}
    dvdy::AbstractArray{T, 2}
    dudy_c::AbstractArray{T, 2}
    dvdx_c::AbstractArray{T, 2}
    shear_c::AbstractArray{T, 2}
    shear_h::AbstractArray{T, 2}
    β_crop::AbstractArray{T, 2}
    gf_crop::AbstractArray{T, 2}
    denu::AbstractArray{T, 2}
    denv::AbstractArray{T, 2}
    ipolgfu::AbstractArray{T, 2}
    ipolgfv::AbstractArray{T, 2}
    hη::AbstractArray{T, 2}
    hη_c::AbstractArray{T, 2}
    rhs::AbstractVector{T}
    f1::AbstractVector{T}
    f2::AbstractVector{T}
    f3::AbstractVector{T}
    sui::AbstractVector{T}
    hui::AbstractVector{T}
    dui::AbstractVector{T}
    svi::AbstractVector{T}
    hvi::AbstractVector{T}
    dvi::AbstractVector{T}
    uvfixed::AbstractVector{T}
    start_guess::AbstractVector{T}
    picard_resid::AbstractVector{T}
    picard_correction::AbstractVector{T}
    gs_resid::AbstractVector{T}
    prolonged::AbstractVector{T}
    op_coarse_tmp1::AbstractVector{T}
    op_coarse_tmp2::AbstractVector{T}
    op_diag::AbstractVector{T}
    r_xx::AbstractArray{T, 2}
    r_yy::AbstractArray{T, 2}
    r_xy::AbstractArray{T, 2}
    extra::AbstractArray{T, 2}
    dx_inv::T
    dy_inv::T
    gs_colour_indices::Vector{AbstractVector{Int}} = AbstractVector{Int}[Int[], Int[], Int[], Int[]]
    gs_colours_filled::Bool = false
    mg_ops::MultigridScratch{T}
    gu_inner_index_map::AbstractArray{Int, 2}
    gv_inner_index_map::AbstractArray{Int, 2}
    haar_u::AbstractArray{T, 2}
    haar_u_tmp::AbstractArray{T, 2}
    haar_v::AbstractArray{T, 2}
    haar_v_tmp::AbstractArray{T, 2}
    glen_a_ref::AbstractArray{T, 2}
end

"""
    stencil_scratch!(model)

Return persistent scratch for the momentum operator and Picard stencil applies.
Allocated on first use. Rebuilt if the cached buffer is on a different
backend to `gh.h`, which happens after a GPU checkpoint pickup.
"""
function stencil_scratch!(model::AbstractModel{T, N}) where {T, N}
    ref = model.fields.stencil_scratch
    s = ref[]
    proto = model.fields.gh.h
    if s !== nothing && KA.get_backend(s.rhs) === KA.get_backend(proto)
        return s::StencilScratch{T}
    end
    allocated = allocate_stencil_scratch(model)
    ref[] = allocated
    return allocated
end

zeros_like(prototype::AbstractArray{T}, dims::Integer...) where {T} =
    fill!(similar(prototype, T, dims...), zero(T))

function copy_like(prototype::AbstractArray, x::AbstractArray)
    y = similar(prototype, eltype(x), size(x)...)
    copyto!(y, x)
    return y
end

"""
    copy_onto!(dest, src)

Write `src` into `dest`. `src` may be a scalar or a host array.

Do not use `dest .= host_matrix` when `dest` is a GPU array: that broadcast
tries to launch a kernel with a `Matrix`, which is not allowed.
"""
copy_onto!(dest::AbstractArray, src::Number) = (fill!(dest, src); dest)
copy_onto!(dest::AbstractArray, src::AbstractArray) = (copyto!(dest, src); dest)

"""
    _host(x)

Copy `x` to a host `Array` if needed. `_host(::Array)` is a no-copy so the
CPU path does not duplicate every field.
"""
_host(x) = Array(x)
_host(x::Array) = x

function allocate_stencil_scratch(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    grid = model.grid
    proto = gh.h

    # Index maps are built on the host (findall / column-major numbering),
    # then copied onto the same backend as `gh.h`.
    mask_u = Array(gu.mask_inner)
    mask_v = Array(gv.mask_inner)
    gu_inner_indices = copy_like(vec(proto), findall(vec(mask_u)))
    gv_inner_indices = copy_like(vec(proto), findall(vec(mask_v)))
    gu_inner_index_map = copy_like(proto, inner_index_map(mask_u))
    gv_inner_index_map = copy_like(proto, inner_index_map(mask_v))

    # Preallocate intermediate variables used by op_fun and Picard applies
    dudx = similar(gh.h)
    dvdy = similar(gh.h)
    r_xx = similar(gh.h)
    r_yy = similar(gh.h)

    dudy_c = similar(gh.h, T, gc.nxc, gc.nyc)
    dvdx_c = similar(dudy_c)
    r_xy = similar(dudy_c)

    extra = similar(gh.h)

    surf_crop = similar(gh.h)
    ones_crop = similar(gh.h)
    tmpu = similar(gu.u)
    tmpv = similar(gv.v)
    tmpui = similar(gu.u, gu.ni)
    tmpvi = similar(gv.v, gv.ni)
    u_crop = similar(gu.u)
    v_crop = similar(gv.v)
    u_h = similar(gh.h)
    v_h = similar(gh.h)
    shear_c = similar(dudy_c)
    shear_h = similar(gh.h)
    β_crop = similar(gh.h)
    gf_crop = similar(gh.h)
    denu = similar(gu.u)
    denv = similar(gv.v)
    ipolgfu = zeros_like(gu.u, gu.nxu, gu.nyu)
    ipolgfv = zeros_like(gv.v, gv.nxv, gv.nyv)
    hη = similar(gh.h)
    hη_c = similar(gh.h, T, gc.nxc, gc.nyc)
    rhs = zeros_like(vec(proto), gu.ni + gv.ni)
    f1 = zeros_like(vec(proto), gu.ni + gv.ni)
    f2 = zeros_like(vec(proto), gu.ni + gv.ni)
    f3 = zeros_like(vec(proto), gu.ni + gv.ni)
    sui = zeros_like(vec(proto), gu.ni)
    hui = zeros_like(vec(proto), gu.ni)
    dui = zeros_like(vec(proto), gu.ni)
    svi = zeros_like(vec(proto), gv.ni)
    hvi = zeros_like(vec(proto), gv.ni)
    dvi = zeros_like(vec(proto), gv.ni)
    uvfixed = zeros_like(vec(proto), gu.nxu * gu.nyu + gv.nxv * gv.nyv)

    # Inner velocity length (free u-points, then free v-points).
    ni = gu.ni + gv.ni

    # Picard / smoother vectors. Reused every iterate instead of similar/zero.
    start_guess = zeros_like(vec(proto), ni)
    picard_resid = zeros_like(vec(proto), ni)
    picard_correction = zeros_like(vec(proto), ni)
    gs_resid = zeros_like(vec(proto), ni)
    prolonged = zeros_like(vec(proto), ni)
    op_coarse_tmp1 = zeros_like(vec(proto), ni)
    op_coarse_tmp2 = zeros_like(vec(proto), ni)

    # Jacobi diagonal (filled each get_op_diag from the fused stencil).
    op_diag = zeros_like(vec(proto), ni)

    dx_inv = one(T) / grid.dx
    dy_inv = one(T) / grid.dy

    # Two work arrays per grid: each Haar step reads one and writes the other.
    haar_u = similar(gu.u)
    haar_u_tmp = similar(gu.u)
    haar_v = similar(gv.v)
    haar_v_tmp = similar(gv.v)
    glen_a_ref = copy_like(proto, model.params.glen_a_ref)

    return StencilScratch{T}(;
        gu_inner_indices,
        gv_inner_indices,
        surf_crop,
        ones_crop,
        tmpu,
        tmpv,
        tmpui,
        tmpvi,
        u_crop,
        v_crop,
        u_h,
        v_h,
        dudx,
        dvdy,
        dudy_c,
        dvdx_c,
        shear_c,
        shear_h,
        β_crop,
        gf_crop,
        denu,
        denv,
        ipolgfu,
        ipolgfv,
        hη,
        hη_c,
        rhs,
        f1,
        f2,
        f3,
        sui,
        hui,
        dui,
        svi,
        hvi,
        dvi,
        uvfixed,
        start_guess,
        picard_resid,
        picard_correction,
        gs_resid,
        prolonged,
        op_coarse_tmp1,
        op_coarse_tmp2,
        op_diag,
        r_xx,
        r_yy,
        r_xy,
        extra,
        dx_inv,
        dy_inv,
        mg_ops = MultigridScratch{T}(;
            n_wu = 0,
            n_wv = 0,
            b_coarse = T[],
            correction_coarse = T[],
        ),
        gu_inner_index_map,
        gv_inner_index_map,
        haar_u,
        haar_u_tmp,
        haar_v,
        haar_v_tmp,
        glen_a_ref,
    )
end

"""
    apply_momentum_op!(out, in, s, gh, gu, gv, gc; vecSampled=true)

Multiply the stacked inner (or full-grid) velocity vector by the momentum operator.
Work arrays live on `s`. Rheology diagonals are read from the grids so Picard
updates are seen without rebuilding scratch.
"""
function apply_momentum_op!(
    opvecprod::AbstractVector,
    inputVector::AbstractVector,
    s::StencilScratch{T},
    gh, gu, gv, gc;
    vecSampled::Bool = true,
) where {T}
    if vecSampled
        @assert length(inputVector) == (gu.ni + gv.ni)
        u_in = view(inputVector, 1:gu.ni)
        v_in = view(inputVector, (gu.ni + 1):(gu.ni + gv.ni))
    else
        nu = gu.nxu * gu.nyu
        nv = gv.nxv * gv.nyv
        @assert length(inputVector) == (nu + nv)
        u_in = view(inputVector, 1:nu)
        v_in = view(inputVector, (nu + 1):(nu + nv))
    end
    out_u = view(opvecprod, 1:gu.ni)
    out_v = view(opvecprod, (gu.ni + 1):(gu.ni + gv.ni))

    D_h = reshape(gh.dneghηav[].diag, gh.nxh, gh.nyh)
    D_c = reshape(gc.dneghηav[].diag, gc.nxc, gc.nyc)
    D_u = reshape(gu.dnegβeff[].diag, gu.nxu, gu.nyu)
    D_v = reshape(gv.dnegβeff[].diag, gv.nxv, gv.nyv)
    D_imp = reshape(gh.dimplicit[].diag, gh.nxh, gh.nyh)

    launch!(
        _op_h_stresses!, s.r_xx, s.r_yy, s.extra, s.r_xy, u_in, v_in,
        s.gu_inner_index_map, s.gv_inner_index_map,
        gu.h, gv.h, gu.mask, gv.mask, gc.mask, D_h, D_imp, D_c, s.dx_inv, s.dy_inv, vecSampled;
        ndrange = size(s.r_xx),
    )
    launch!(
        _op_force_u!, out_u, s.r_xx, s.r_xy, s.extra, u_in, s.gu_inner_index_map, gu.h,
        D_u, s.dx_inv, s.dy_inv, vecSampled;
        ndrange = size(s.gu_inner_index_map), sync = false,
    )
    launch!(
        _op_force_v!, out_v, s.r_yy, s.r_xy, s.extra, v_in, s.gv_inner_index_map, gv.h,
        D_v, s.dx_inv, s.dy_inv, vecSampled;
        ndrange = size(s.gv_inner_index_map), sync = false,
    )
    KA.synchronize(KA.get_backend(s.r_xx))
    return opvecprod
end

"""
    get_op_fun(model::AbstractModel)

Returns a function that multiplies a vector by the momentum operator.
"""
function get_op_fun(model::AbstractModel)
    s = stencil_scratch!(model)
    gh, gu, gv, gc = model.fields.gh, model.fields.gu, model.fields.gv, model.fields.gc
    return (out, in; vecSampled = true) -> apply_momentum_op!(out, in, s, gh, gu, gv, gc; vecSampled)
end

"""
    haar_steps(levels)

Return the pairing gaps for a Haar transform with the given number of levels.
The first gap is 2 cells, then 4, then 8, doubling until `levels + 1` gaps
have been listed.
"""
haar_steps(levels::Integer) = ntuple(i -> 2^i, levels + 1)

# KernelAbstractions schedules CPU threads one workgroup at a time. A
# default workgroup that covers the whole ndrange runs on one thread.
# On the main thread, when Julia has extra threads, use one item
# per workgroup so the groups can map onto those threads. Serial, MPI,
# and nested ThreadedSpec workers keep a single workgroup.
function _haar_workgroup()
    (Threads.nthreads() > 1 && Threads.threadid() == 1) ? 1 : nothing
end

"""
    haar_idwt!(a, b, levels, transpose=false) -> result

Apply the Haar inverse transform to the grid `a` by pairing neighbours at
spacings 2, 4, 8, ... instead of a matrix product.
`b` is a workspace of the same size. Each axis reads one array and writes the
other. Pairings run along y first, then along x. The array that contains the
result is returned. If `transpose` is true, this is the adjoint used by
restrict (see `haar_idwtᵀ!`).
"""
function haar_idwt!(a, b, levels, transpose=false)
    steps = haar_steps(levels)
    step_iter = transpose ? steps : reverse(steps)
    return _haar_lift_axes!(a, b, step_iter, transpose)
end

"""
    haar_dwt!(a, b, levels) -> result

Apply the forward Haar transform to the grid `a`.
`update_wavelets!` uses this to form wavelet coefficients from velocity.
`b` is a workspace of the same size. Matches `wavelet_matrix(..., "forward")`.
This is not the reverse Haar used by prolong (`haar_idwt!`).
"""
function haar_dwt!(a, b, levels)
    return _haar_lift_axes!(a, b, haar_steps(levels), Val{:forward}())
end

function _haar_lift_axes!(a, b, step_iter, mix)
    if KA.get_backend(a) isa KA.CPU
        return _haar_lift_axes_line!(a, b, step_iter, mix)
    end
    return _haar_lift_axes_2d!(a, b, step_iter, mix)
end

"""
    _haar_lift_axes_line!(a, b, step_iter, mix)

Apply every Haar spacing along y, then along x, with one fused line kernel
per axis. Used on CPU so RAP does not launch a full-grid kernel per level.
"""
function _haar_lift_axes_line!(a, b, step_iter, mix)
    src, dst = a, b
    nx = size(a, 1)
    wg = _haar_workgroup()
    nwork = wg === 1 ? min(nx, Threads.nthreads()) : 1
    launch!(_haar_lift_y_all!, src, dst, step_iter, mix, nwork;
            ndrange = nwork, workgroupsize = wg)
    if isodd(length(step_iter))
        src, dst = dst, src
    end
    launch!(_haar_lift_x_all!, src, dst, step_iter, mix;
            ndrange = size(a, 2), workgroupsize = wg)
    if isodd(length(step_iter))
        src, dst = dst, src
    end
    return src
end

"""
    _haar_lift_axes_2d!(a, b, step_iter, mix)

Apply each Haar spacing as a 2D kernel (`ndrange = (nx, ny)`).
Used on GPU so each pairing fills the device. CPU keeps `_haar_lift_axes_line!`.
Queue every spacing, then synchronise once so the host does not wait per level.
"""
function _haar_lift_axes_2d!(a, b, step_iter, mix)
    src, dst = a, b
    ndrange = size(a)
    backend = KA.get_backend(a)
    for step in step_iter
        launch!(_haar_lift_y!, dst, src, step, mix; ndrange = ndrange, sync = false)
        src, dst = dst, src
    end
    for step in step_iter
        launch!(_haar_lift_x!, dst, src, step, mix; ndrange = ndrange, sync = false)
        src, dst = dst, src
    end
    KA.synchronize(backend)
    return src
end

"""
    haar_idwtᵀ!(a, b, levels) -> result

Apply the adjoint Haar transform.
Restrict uses this to form wavelet coefficients.
"""
haar_idwtᵀ!(a, b, levels) = haar_idwt!(a, b, levels, true)

"""
    get_restrict_fun(model::AbstractModel)

Map a residual on free ice points onto the coarse wavelet coefficients
kept above the threshold. Matches the old sparse `samp * idwtᵀ * spread_inner`.
"""
function get_restrict_fun(model::AbstractModel)
    s = stencil_scratch!(model)
    @unpack wu, wv, gu, gv = model.fields

    function restrict_fun!(restrictvec::AbstractVector, vec::AbstractVector)
        @assert length(vec) == (gu.ni + gv.ni)
        n_wu = wu.n[]
        n_wv = wv.n[]

        # spread_inner
        launch!(_scatter_mapped!, s.haar_u, view(vec, 1:gu.ni), s.gu_inner_index_map;
                ndrange = size(s.haar_u), sync = false)
        launch!(_scatter_mapped!, s.haar_v, view(vec, (gu.ni + 1):(gu.ni + gv.ni)), s.gv_inner_index_map;
                ndrange = size(s.haar_v))

        # idwtᵀ
        ru = haar_idwtᵀ!(s.haar_u, s.haar_u_tmp, wu.levels)
        rv = haar_idwtᵀ!(s.haar_v, s.haar_v_tmp, wv.levels)

        # samp
        launch!(_gather_mapped!, view(restrictvec, 1:n_wu), ru, wu.index_map;
                ndrange = size(wu.index_map), sync = false)
        launch!(_gather_mapped!, view(restrictvec, (n_wu + 1):(n_wu + n_wv)), rv, wv.index_map;
                ndrange = size(wv.index_map))
        return restrictvec
    end

    return restrict_fun!
end

"""
    get_prolong_fun(model::AbstractModel)

Map kept coarse wavelet coefficients back to a velocity increment on free
ice points. Matches the old sparse `samp_inner * idwt * spread`.
"""
function get_prolong_fun(model::AbstractModel)
    s = stencil_scratch!(model)
    @unpack wu, wv, gu, gv = model.fields

    function prolong_fun!(prolongvec::AbstractVector, waveletvec::AbstractVector)
        n_wu = wu.n[]
        n_wv = wv.n[]
        @assert length(waveletvec) == (n_wu + n_wv)

        # spread
        launch!(_scatter_mapped!, s.haar_u, view(waveletvec, 1:n_wu), wu.index_map;
                ndrange = size(wu.index_map), sync = false)
        launch!(_scatter_mapped!, s.haar_v, view(waveletvec, (n_wu + 1):(n_wu + n_wv)), wv.index_map;
                ndrange = size(wv.index_map))

        # idwt
        ru = haar_idwt!(s.haar_u, s.haar_u_tmp, wu.levels)
        rv = haar_idwt!(s.haar_v, s.haar_v_tmp, wv.levels)

        # samp_inner
        launch!(_gather_mapped!, view(prolongvec, 1:gu.ni), ru, s.gu_inner_index_map;
                ndrange = size(s.gu_inner_index_map), sync = false)
        launch!(_gather_mapped!, view(prolongvec, (gu.ni + 1):(gu.ni + gv.ni)), rv, s.gv_inner_index_map;
                ndrange = size(s.gv_inner_index_map))
        return prolongvec
    end

    return prolong_fun!
end

"""
    _quadrant_positive_area(z1, z2, z3)

How much of one cell quarter sits above zero (0 to 1).
`z1` is this cell; `z2` and `z3` are the two neighbours for this quarter.
"""
@inline function _quadrant_positive_area(z1::T, z2::T, z3::T) where {T}
    if (z1 > 0) && (z2 > 0) && (z3 > 0)
        return one(T)
    end
    if (z1 <= 0) && (z2 <= 0) && (z3 <= 0)
        return zero(T)
    end

    # Pick the neighbour with the larger change so we do not divide by almost zero.
    flip = abs(z3 - z1) < abs(z2 - z1)
    den = flip ? (z2 - z1) : (z3 - z1)
    a = -(flip ? (z3 - z1) : (z2 - z1)) / den
    b = -T(2) * z1 / den
    a1 = T(0.5) * (b * b) / a
    a2 = T(0.5) * ((one(T) - b) * (one(T) - b)) / a
    a3 = T(0.5) * a + b
    test1 = Int(b > 0) + Int(b > 1)
    test2 = Int((a + b) > 0) + Int((a + b) > 1)
    ix = 1 + test1 + 3 * test2
    area = if ix == 1
        zero(T)
    elseif ix == 2
        -a1
    elseif ix == 3
        a2 - a1
    elseif ix == 4
        a1 + a3
    elseif ix == 5
        a3
    elseif ix == 6
        a2 + a3
    elseif ix == 7
        one(T) - a2 + a1
    elseif ix == 8
        one(T) - a2
    else
        one(T)
    end
    if b < 0
        area = one(T) - area
    end
    if z1 < 0
        area = one(T) - area
    end
    return area
end

@inline function _neighbour_z(z, mask, i, j, i_nb, j_nb)
    @inbounds mask[i_nb, j_nb] ? z[i_nb, j_nb] : z[i, j]
end

"""
    _quadrant_area(z, mask, i, j, m, n, quadrant)

Above-zero area of one quarter of cell `(i, j)`.
Quadrant 1 is east and north, 2 north and west, 3 west and south, 4 south and east.
Invalid cells return 0. Invalid or off-grid neighbours reuse this cell's value.
"""
@inline function _quadrant_area(z, mask, i, j, m, n, quadrant)
    @inbounds begin
        mask[i, j] || return zero(eltype(z))
        z1 = z[i, j]
        if quadrant == 1
            z2 = _neighbour_z(z, mask, i, j, min(i + 1, m), j)
            z3 = _neighbour_z(z, mask, i, j, i, min(j + 1, n))
        elseif quadrant == 2
            z2 = _neighbour_z(z, mask, i, j, i, min(j + 1, n))
            z3 = _neighbour_z(z, mask, i, j, max(i - 1, 1), j)
        elseif quadrant == 3
            z2 = _neighbour_z(z, mask, i, j, max(i - 1, 1), j)
            z3 = _neighbour_z(z, mask, i, j, i, max(j - 1, 1))
        else
            z2 = _neighbour_z(z, mask, i, j, i, max(j - 1, 1))
            z3 = _neighbour_z(z, mask, i, j, min(i + 1, m), j)
        end
        return _quadrant_positive_area(z1, z2, z3)
    end
end

@inline function _cell_quadrants(z, mask, i, j, m, n)
    return (
        _quadrant_area(z, mask, i, j, m, n, 1),
        _quadrant_area(z, mask, i, j, m, n, 2),
        _quadrant_area(z, mask, i, j, m, n, 3),
        _quadrant_area(z, mask, i, j, m, n, 4),
    )
end

# This cell's four quarters, then u from the west neighbour and v from the south.
@kernel function _pos_fraction_kernel!(area_h, area_u, area_v, z, mask)
    i, j = @index(Global, NTuple)
    m, n = size(z)
    T = eltype(area_h)
    quarter = T(0.25)
    @inbounds begin
        q1 = zero(T)
        q2 = zero(T)
        q3 = zero(T)
        q4 = zero(T)
        if (i <= m) && (j <= n)
            q1, q2, q3, q4 = _cell_quadrants(z, mask, i, j, m, n)
            area_h[i, j] = quarter * (q1 + q2 + q3 + q4)
        end

        if j <= n
            acc = (i <= m) ? quarter * (q2 + q3) : zero(T)
            if i > 1
                w1, w2, w3, w4 = _cell_quadrants(z, mask, i - 1, j, m, n)
                acc += quarter * (w1 + w4)
            end
            area_u[i, j] = acc
        end

        if i <= m
            acc = (j <= n) ? quarter * (q3 + q4) : zero(T)
            if j > 1
                s1, s2, s3, s4 = _cell_quadrants(z, mask, i, j - 1, m, n)
                acc += quarter * (s1 + s2)
            end
            area_v[i, j] = acc
        end
    end
end

function _pos_fraction_mask(z, mask)
    size(mask) == size(z) || throw(DimensionMismatch(
        "pos_fraction mask must have size $(size(z)), got $(size(mask))"
    ))
    host_mask = mask isa BitArray ? Array(mask) : mask
    if typeof(host_mask) == typeof(similar(z, Bool))
        return host_mask
    end
    msk = similar(z, Bool)
    copyto!(msk, host_mask)
    return msk
end

"""
    pos_fraction!(area_h, area_u, area_v, z1, mask)

In-place `pos_fraction`. Inputs must all be on the CPU, or all on the GPU.
"""
function pos_fraction!(area_h, area_u, area_v, z1, mask)
    m, n = size(z1)
    size(area_h) == (m, n) || throw(DimensionMismatch("area_h must be $((m, n))"))
    size(area_u) == (m + 1, n) || throw(DimensionMismatch("area_u must be $((m + 1, n))"))
    size(area_v) == (m, n + 1) || throw(DimensionMismatch("area_v must be $((m, n + 1))"))
    msk = _pos_fraction_mask(z1, mask)
    launch!(_pos_fraction_kernel!, area_h, area_u, area_v, z1, msk; ndrange = (m + 1, n + 1))
    return area_h, area_u, area_v
end

"""
pos_fraction(z1;mask=mask) -> area_fraction, area_fraction_u, area_fraction_v

Return fraction of each grid cell with function z1 above zero. Uses bilinear
interpolation of values at three nearest cell centers to represent the function.
In:
   z1:             m x n array of gridded function values.
   mask:           m x n mask   1 = valid data, 0= invalid data.
Out:
   area_fraction:   m x n array showing area fraction of interpolated z1>0 on h-grid.
   area_fraction_u: m+1 x n array showing area fraction of interpolated z1>0 on u-grid.
   area_fraction_v: m x n+1 array showing area fraction of interpolated z1>0 on v-grid.
"""
function pos_fraction(z1::AbstractArray{T,2}; mask=trues(size(z1))) where {T}
    m, n = size(z1)
    area_h = similar(z1)
    area_u = similar(z1, T, m + 1, n)
    area_v = similar(z1, T, m, n + 1)
    return pos_fraction!(area_h, area_u, area_v, z1, mask)
end



#MISMIP+ bed elevation
function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75
    wc = 24000.0; fc = 4000.0; dc = 500.0
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) )
    b = max(bx(x) + by(y), -720.0)
    return b
end

"""
     get_glx(model)

Return the grounding line in the form x = x(y). Assumes each y-row has at least
one grid point where ice grounded and one where ice floating.
"""
function get_glx(model)
    @unpack fields, grid = model
    glmask=diff(sign.(fields.gh.haf),dims=1).==-2 #calculate where sign of height above floating passes thru zero
    glx1=grid.xxh[1:end-1,:][glmask] #x co-ordiates upstream of grounding line
    glx2=grid.xxh[2:end,:][glmask] #x co-ordinates immediately downstream
    haf1=fields.gh.haf[1:end-1,:][glmask] #Height above floatation immediately upstream of gl
    haf2=fields.gh.haf[2:end,:][glmask]
    glx=glx1+(glx2-glx1).*(zero(haf1)-haf1)./(haf2-haf1) #interpolate between grid points to find grounding line
    return glx
end




"""
    icedraft(s,h,sea_level_wrt_geoid)

Compute ice draft.
"""
icedraft(s,h,sea_level_wrt_geoid)=max(h-(s-sea_level_wrt_geoid),zero(typeof(h-(s-sea_level_wrt_geoid))))


"""
    height_above_floatation(h,b,params)

Compute height above floatation.

"""
height_above_floatation(h,b,params) = h - (params.density_ocean/params.density_ice)*(params.sea_level_wrt_geoid - b)


"""
    volume_above_floatation(h,b,params)

Compute the volume above floatation: integrated height above floatation for cells with positive height above floatation
"""
volume_above_floatation(h,b,params,grid) = sum(sum(height_above_floatation.(h,b,params)[height_above_floatation.(h,b, params) .> 0])) .* grid.dx .* grid.dy

"""
    glen_b(temperature,damage,params)

Compute stiffness parameter B in Glen flow law.

"""

function glen_b(temperature,damage,glen_a_ref, glen_n, glen_a_activation_energy, glen_temperature_ref, gas_const)
    glen_a0 = glen_a_ref*exp(+glen_a_activation_energy/(glen_temperature_ref*gas_const) )
    glen_b = (1-damage)*( glen_a0*exp(-glen_a_activation_energy/(temperature*gas_const)) )^(-1.0/glen_n)
    return glen_b
end

@kernel function _update_glen_b_kernel!(
    glen_b_arr,
    θ,
    Φ,
    glen_a_ref,
    glen_n,
    glen_a_activation_energy,
    glen_temperature_ref,
    gas_const,
)
    i, j, k = @index(Global, NTuple)
    @inbounds glen_b_arr[i, j, k] = glen_b(
        θ[i, j, k],
        Φ[i, j, k],
        glen_a_ref[i, j],
        glen_n,
        glen_a_activation_energy,
        glen_temperature_ref,
        gas_const,
    )
end

"""
    fill_glen_b!(glen_b_arr, θ, Φ, glen_a_ref, glen_n, glen_a_activation_energy, glen_temperature_ref, gas_const)

Fill `glen_b_arr` from temperature and damage using Glen's flow-law formula.

This is the same calculation as `update_glen_b!`. Model construction calls it so the
starting Glen B field matches what later timesteps compute.
"""
function fill_glen_b!(
    glen_b_arr,
    θ,
    Φ,
    glen_a_ref,
    glen_n,
    glen_a_activation_energy,
    glen_temperature_ref,
    gas_const,
)
    a_ref = glen_a_ref
    if !_is_host_array(glen_b_arr) && _is_host_array(glen_a_ref)
        a_ref = similar(glen_b_arr, eltype(glen_a_ref), size(glen_a_ref))
        copyto!(a_ref, glen_a_ref)
    end
    launch!(
        _update_glen_b_kernel!,
        glen_b_arr,
        θ,
        Φ,
        a_ref,
        glen_n,
        glen_a_activation_energy,
        glen_temperature_ref,
        gas_const;
        ndrange = size(glen_b_arr),
    )
    return glen_b_arr
end


"""
    get_u_mask(h_mask)

Find mask of valid grid points on u-grid corresponding to a mask defined on h-grid.

"""
function get_u_mask(h_mask)
    #include all u faces next to a selected center
    (nx,ny)=size(h_mask)
    u_mask=falses(nx+1,ny)
    u_mask[1:end-1,1:end]=u_mask[1:end-1,1:end].|h_mask
    u_mask[2:end,1:end]=u_mask[2:end,1:end].|h_mask
    return u_mask
end
"""
    get_v_mask(h_mask)

Find mask of valid grid points on v-grid corresponding to a mask defined on h-grid.

"""
function get_v_mask(h_mask)
    #include all v faces next to a selected center
    (nx,ny)=size(h_mask)
    v_mask=falses(nx,ny+1)
    v_mask[1:end,1:end-1]=v_mask[1:end,1:end-1].|h_mask
    v_mask[1:end,2:end]=v_mask[1:end,2:end].|h_mask
    return v_mask
end
"""
    get_c_mask(h_mask)

Find mask of valid grid points on c-grid corresponding to a mask defined on h-grid.

"""
function get_c_mask(h_mask)
    #select cell corners with four neighbouring cell centers in h_mask
    c_mask=h_mask[1:end-1,1:end-1] .& h_mask[1:end-1,2:end] .& h_mask[2:end,1:end-1] .& h_mask[2:end,2:end]
    return c_mask
end
"""
    clip(trial_mask)

Find mask of valid grid points on h-grid corresponding to a trial mask, also defined on h-grid.
Clip any grid points from the trial mask that cannot be used in the model.

"""
function clip(trial_mask)

    good_corners=get_c_mask(trial_mask)

    #include all centers next to a selected corner
    (nx,ny)=size(trial_mask)
    mask=falses(nx,ny)
    mask[1:end-1,1:end-1]=mask[1:end-1,1:end-1].|good_corners
    mask[1:end-1,2:end]=mask[1:end-1,2:end].|good_corners
    mask[2:end,1:end-1]=mask[2:end,1:end-1].|good_corners
    mask[2:end,2:end]=mask[2:end,2:end].|good_corners

    return mask
end


"""
     get_resid(x,op,b)

Function to return residual b - op x

"""
function get_resid(x,op,b)
    resid=similar(b)
    get_resid!(resid,x,op,b)
end


"""
     get_resid!(resid,x,op,b)

In-place function to return residual b - op x

"""
function get_resid!(resid,x,op,b)
    mem_resid= @allocated   mul!(resid,op,x)
   # println("Memory allocated in get_resid is: ", mem_resid, " bytes")
    resid .= b .- resid
end

end
