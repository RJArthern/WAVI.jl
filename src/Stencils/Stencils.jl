module Stencils

export _diff_x!,
    _diff_y!,
    _diff_xT!,
    _diff_yT!,
    _diff_x_staggered!,
    _diff_y_staggered!,
    _diff_xT_staggered!,
    _diff_yT_staggered!,
    _avg_x!,
    _avg_y!,
    _avg_xT!,
    _avg_yT!,
    _avg_xy!,
    _avg_xyT!,
    _apply_mask!,
    _scale!,
    _scale_sum!,
    _add_scale!,
    _masked_mul!,
    _masked_scale_sum!,
    _gather!,
    _scatter!,
    _scatter_mapped!,
    _gather_mapped!,
    _haar_lift_x!,
    _haar_lift_y!,
    _op_h_stresses!,
    _op_force_u!,
    _op_force_v!,
    _op_diag_u!,
    _op_diag_v!,
    launch!

using KernelAbstractions: KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

"""
    launch!(kernel!, args...; ndrange, sync = true, workgroupsize = nothing)

Wrapper to launch a KernelAbstractions `kernel!`.
The backend is taken from the first kernel argument via `KA.get_backend`.
On a CPU with more than one Julia thread, `CPU(static = true)` is used so each
kernel does not pay an `@spawn` per launch. Nested calls from other threads
keep the default dynamic CPU backend.
Ensures synchronisation after the kernel is launched unless `sync = false`.
"""
function launch!(kernel!, args...; ndrange, sync = true, workgroupsize = nothing)
    backend = KA.get_backend(first(args))
    # Static assignment on the main thread only. Nested ThreadedSpec
    # workers keep dynamic spawn so they do not nest `@threads :static`.
    if backend isa KA.CPU && Threads.nthreads() > 1 && Threads.threadid() == 1
        backend = KA.CPU(static = true)
    end
    if workgroupsize === nothing
        kernel!(backend)(args...; ndrange = ndrange)
    else
        kernel!(backend)(args...; ndrange = ndrange, workgroupsize = workgroupsize)
    end
    if sync
        KA.synchronize(backend)
    end
end

# Finite Differences

"""
    _diff_x!(out, inp, dx_inv)

Compute forward finite difference of `inp` in x-direction.
Multiplies the result by `dx_inv` (which is `1/dx`) and stores in `out`.
Equivalent to `out = ∂x * inp` in sparse matrix notation.
"""
@kernel function _diff_x!(out, inp, dx_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i + 1, j] - inp[i, j]) * dx_inv
end

"""
    _diff_y!(out, inp, dy_inv)

Compute forward finite difference of `inp` in y-direction.
Multiplies the result by `dy_inv` (which is `1/dy`) and stores in `out`.
Equivalent to `out = ∂y * inp` in sparse matrix notation.
"""
@kernel function _diff_y!(out, inp, dy_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i, j + 1] - inp[i, j]) * dy_inv
end

"""
    _diff_xT!(out, inp, dx_inv)

Compute backward (transpose) finite difference of `inp` in x-direction.
Equivalent to `out = ∂xᵀ * inp` in sparse matrix notation.
`out` is one cell wider than `inp` (U-grid); edge values are treated as zero.
"""
@kernel function _diff_xT!(out, inp, dx_inv)
    i, j = @index(Global, NTuple)
    nx = size(inp, 1)
    val_left = i > 1 ? inp[i - 1, j] : zero(eltype(inp))
    val_right = i <= nx ? inp[i, j] : zero(eltype(inp))
    @inbounds out[i, j] = (val_left - val_right) * dx_inv
end

"""
    _diff_yT!(out, inp, dy_inv)

Compute backward (transpose) finite difference of `inp` in y-direction.
Equivalent to `out = ∂yᵀ * inp` in sparse matrix notation.
`out` is one cell taller than `inp` (V-grid); edge values are treated as zero.
"""
@kernel function _diff_yT!(out, inp, dy_inv)
    i, j = @index(Global, NTuple)
    ny = size(inp, 2)
    val_bot = j > 1 ? inp[i, j - 1] : zero(eltype(inp))
    val_top = j <= ny ? inp[i, j] : zero(eltype(inp))
    @inbounds out[i, j] = (val_bot - val_top) * dy_inv
end

# Staggered Finite Differences (C-grid shear)

"""
    _diff_y_staggered!(out, inp, dy_inv)

Forward difference of U-grid `inp` in y onto the C-grid, using interior faces.
Equivalent to `out = gu.∂y * inp` (`∂1d ⊗ χ`) in sparse matrix notation.
"""
@kernel function _diff_y_staggered!(out, inp, dy_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i + 1, j + 1] - inp[i + 1, j]) * dy_inv
end

"""
    _diff_x_staggered!(out, inp, dx_inv)

Forward difference of V-grid `inp` in x onto the C-grid, using interior faces.
Equivalent to `out = gv.∂x * inp` (`χ ⊗ ∂1d`) in sparse matrix notation.
"""
@kernel function _diff_x_staggered!(out, inp, dx_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i + 1, j + 1] - inp[i, j + 1]) * dx_inv
end

"""
    _diff_yT_staggered!(out, inp, dy_inv)

Transpose of `_diff_y_staggered!`. Equivalent to `out = gu.∂yᵀ * inp`.
"""
@kernel function _diff_yT_staggered!(out, inp, dy_inv)
    i, j = @index(Global, NTuple)
    nxc = size(inp, 1)
    nyc = size(inp, 2)
    if i > 1 && i <= nxc + 1
        val_bot = j > 1 ? inp[i - 1, j - 1] : zero(eltype(inp))
        val_top = j <= nyc ? inp[i - 1, j] : zero(eltype(inp))
        @inbounds out[i, j] = (val_bot - val_top) * dy_inv
    else
        @inbounds out[i, j] = zero(eltype(out))
    end
end

"""
    _diff_xT_staggered!(out, inp, dx_inv)

Transpose of `_diff_x_staggered!`. Equivalent to `out = gv.∂xᵀ * inp`.
"""
@kernel function _diff_xT_staggered!(out, inp, dx_inv)
    i, j = @index(Global, NTuple)
    nxc = size(inp, 1)
    nyc = size(inp, 2)
    if j > 1 && j <= nyc + 1
        val_left = i > 1 ? inp[i - 1, j - 1] : zero(eltype(inp))
        val_right = i <= nxc ? inp[i, j - 1] : zero(eltype(inp))
        @inbounds out[i, j] = (val_left - val_right) * dx_inv
    else
        @inbounds out[i, j] = zero(eltype(out))
    end
end

# Averaging

"""
    _avg_x!(out, inp)

Compute arithmetic mean of adjacent cells in x-direction.
Equivalent to `out = cent * inp` for u-grid in sparse matrix notation.

Moves a value from cell center to the edge by averaging.
"""
@kernel function _avg_x!(out, inp)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i, j] + inp[i + 1, j]) * 0.5
end

"""
    _avg_y!(out, inp)

Compute arithmetic mean of adjacent cells in y-direction.
Equivalent to `out = cent * inp` for v-grid in sparse matrix notation.

Moves value from cell center to edge by averaging.
"""
@kernel function _avg_y!(out, inp)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i, j] + inp[i, j + 1]) * 0.5
end

"""
    _avg_xT!(out, inp)

Transpose of `_avg_x!`. Equivalent to `out = gu.centᵀ * inp` (H-grid to U-grid).
Edge values are treated as zero.
"""
@kernel function _avg_xT!(out, inp)
    i, j = @index(Global, NTuple)
    nx = size(inp, 1)
    val_left = i > 1 ? inp[i - 1, j] : zero(eltype(inp))
    val_right = i <= nx ? inp[i, j] : zero(eltype(inp))
    @inbounds out[i, j] = (val_left + val_right) * 0.5
end

"""
    _avg_yT!(out, inp)

Transpose of `_avg_y!`. Equivalent to `out = gv.centᵀ * inp` (H-grid to V-grid).
Edge values are treated as zero.
"""
@kernel function _avg_yT!(out, inp)
    i, j = @index(Global, NTuple)
    ny = size(inp, 2)
    val_bot = j > 1 ? inp[i, j - 1] : zero(eltype(inp))
    val_top = j <= ny ? inp[i, j] : zero(eltype(inp))
    @inbounds out[i, j] = (val_bot + val_top) * 0.5
end

"""
    _avg_xy!(out, inp)

Average four neighbouring H-grid cells onto the C-grid.
Equivalent to `out = gh.cent_xy * inp`.
"""
@kernel function _avg_xy!(out, inp)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i, j] + inp[i + 1, j] + inp[i, j + 1] + inp[i + 1, j + 1]) * 0.25
end

"""
    _avg_xyT!(out, inp)

Transpose of `_avg_xy!`. Equivalent to `out = gc.cent * inp` (C-grid to H-grid).
Edge values are treated as zero.
"""
@kernel function _avg_xyT!(out, inp)
    i, j = @index(Global, NTuple)
    nxc = size(inp, 1)
    nyc = size(inp, 2)
    val = zero(eltype(out))
    if i > 1 && j > 1
        val += inp[i - 1, j - 1]
    end
    if i <= nxc && j > 1
        val += inp[i, j - 1]
    end
    if i > 1 && j <= nyc
        val += inp[i - 1, j]
    end
    if i <= nxc && j <= nyc
        val += inp[i, j]
    end
    @inbounds out[i, j] = val * 0.25
end

# Masking & Scaling

"""
    _apply_mask!(arr, mask)

Zero out elements of `arr` where `mask` is false. Leaves elements unchanged
where `mask` is true. Operates in-place on `arr`.
Equivalent to `out = crop * inp`.

Zeroes out any cells that are empty ocean or rock (not ice).
"""
@kernel function _apply_mask!(arr, mask)
    i, j = @index(Global, NTuple)
    @inbounds arr[i, j] = mask[i, j] ? arr[i, j] : zero(eltype(arr))
end

"""
    _scale!(out, inp, diag, factor)

Multiply each element of `inp` by `diag` and `factor`.
Equivalent to `out = factor * Diagonal(diag) * inp`.
"""
@kernel function _scale!(out, inp, diag, factor)
    i = @index(Global, Linear)
    @inbounds out[i] = factor * diag[i] * inp[i]
end

"""
    _scale_sum!(out, x, y, c1, c2, diag, factor)

`out = factor * diag * (c1 * x + c2 * y)`.
"""
@kernel function _scale_sum!(out, x, y, c1, c2, diag, factor)
    i = @index(Global, Linear)
    @inbounds out[i] = factor * diag[i] * (c1 * x[i] + c2 * y[i])
end

"""
    _add_scale!(out, x, y, diag)

`out = diag * (x + y)`.
"""
@kernel function _add_scale!(out, x, y, diag)
    i = @index(Global, Linear)
    @inbounds out[i] = diag[i] * (x[i] + y[i])
end

"""
    _masked_mul!(out, a, b, mask)

`out = mask ? a * b : 0`.
"""
@kernel function _masked_mul!(out, a, b, mask)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = mask[i, j] ? a[i, j] * b[i, j] : zero(eltype(out))
end

"""
    _masked_scale_sum!(out, x, y, diag, mask)

`out = mask ? -diag * (x + y) : 0`.
"""
@kernel function _masked_scale_sum!(out, x, y, diag, mask)
    i = @index(Global, Linear)
    @inbounds out[i] = mask[i] ? -diag[i] * (x[i] + y[i]) : zero(eltype(out))
end

# Scatter / Gather

"""
    _gather!(out_vec, inp_2d, indices)

Extract values from a 2D array `inp_2d` at the linear indices specified by `indices`.
Stores the result compactly in a 1D vector `out_vec`.
Equivalent to `out = samp * inp`.

Takes full 2D map and packs only active ice cells into a small 1D vector for CG solver.
"""
@kernel function _gather!(out_vec, inp_2d, indices)
    k = @index(Global, Linear)
    @inbounds begin
        idx = indices[k]
        out_vec[k] = inp_2d[idx]
    end
end

"""
    _scatter!(out_2d, inp_vec, indices)

Place values from a 1D vector `inp_vec` into a 2D array `out_2d` at the linear
indices specified by `indices`.
Equivalent to `out = spread * inp`.

Takes small 1D CG solver vector and spreads answers back onto full 2D map.
"""
@kernel function _scatter!(out_2d, inp_vec, indices)
    k = @index(Global, Linear)
    @inbounds begin
        idx = indices[k]
        out_2d[idx] = inp_vec[k]
    end
end

"""
    _scatter_mapped!(out_2d, inp_vec, index_map)

Put a short packed vector onto a 2D grid using `index_map`.
Where `index_map[i, j]` is `k > 0`, write `inp_vec[k]`. Where it is `0`, write zero.
Equivalent to `out = spread * inp`, and zeros cells that spread would leave untouched.
"""
@kernel function _scatter_mapped!(out_2d, inp_vec, index_map)
    i, j = @index(Global, NTuple)
    @inbounds begin
        k = index_map[i, j]
        out_2d[i, j] = k == 0 ? zero(eltype(out_2d)) : inp_vec[k]
    end
end

"""
    _gather_mapped!(out_vec, inp_2d, index_map)

Copy a 2D grid into a short packed vector using `index_map`.
Where `index_map[i, j]` is `k > 0`, write `inp_2d[i, j]` into `out_vec[k]`.
Equivalent to `out = samp * inp`.
"""
@kernel function _gather_mapped!(out_vec, inp_2d, index_map)
    i, j = @index(Global, NTuple)
    @inbounds begin
        k = index_map[i, j]
        if k != 0
            out_vec[k] = inp_2d[i, j]
        end
    end
end

# Haar lifting (same even/odd pairing as wavelet_matrix)

"""
    _haar_lift_x!(out, inp, step, transpose)

One Haar pairing along x at spacing `step` (2, 4, 8, ...).
Same even/odd split as `wavelet_matrix(..., "reverse")`.
`idwt` / prolong: even' = even + odd, odd' = odd - even.
`idwtᵀ` / restrict (`transpose=true`) swaps those. Unpaired points are copied.
"""
@kernel function _haar_lift_x!(out, inp, step, transpose)
    i, j = @index(Global, NTuple)
    @inbounds begin
        n = size(inp, 1)
        half = div(step, 2)
        rem = (i - 1) % step
        if rem == 0 && i + half <= n
            odd = inp[i, j]
            even = inp[i + half, j]
            out[i, j] = transpose ? even + odd : odd - even
        elseif rem == half && i > half
            odd = inp[i - half, j]
            even = inp[i, j]
            out[i, j] = transpose ? even - odd : even + odd
        else
            out[i, j] = inp[i, j]
        end
    end
end

"""
    _haar_lift_y!(out, inp, step, transpose)

Same as `_haar_lift_x!`, but along y.
"""
@kernel function _haar_lift_y!(out, inp, step, transpose)
    i, j = @index(Global, NTuple)
    @inbounds begin
        n = size(inp, 2)
        half = div(step, 2)
        rem = (j - 1) % step
        if rem == 0 && j + half <= n
            odd = inp[i, j]
            even = inp[i, j + half]
            out[i, j] = transpose ? even + odd : odd - even
        elseif rem == half && j > half
            odd = inp[i, j - half]
            even = inp[i, j]
            out[i, j] = transpose ? even - odd : even + odd
        else
            out[i, j] = inp[i, j]
        end
    end
end

# Fused momentum operator (SSA resistive stresses and forces)

"""
    velocity_at(u, idx, i, j, vecSampled, z)

Look up the velocity at grid point (i, j).

On the usual solver path, `vecSampled` is true and `u` is the short inner
list. `idx` comes from `inner_index_map`. If this point is inner, the
velocity is `u[idx[i, j]]`. If `idx[i, j]` is `0`, the point is ocean or a
fixed boundary, so return `z` (zero).

On the Dirichlet path, `vecSampled` is false and `u` already holds every
grid point, column by column. Return the entry for (i, j) with no lookup.
"""
@inline function velocity_at(u, idx, i, j, vecSampled, z)
    @inbounds if vecSampled
        k = idx[i, j]
        return k == 0 ? z : u[k]
    else
        return u[i + (j - 1) * size(idx, 1)]
    end
end

"""
    _op_h_stresses!(r_xx, r_yy, extra, r_xy, u, v, u_idx, v_idx, hu, hv, mask_u, mask_v, mask_c, D_h, D_imp, D_c, dx_inv, dy_inv, vecSampled)

On each thickness cell, write the SSA stresses using neighbour velocities
from `velocity_at`.

`r_xx` and `r_yy` are the extensional stresses, and `r_xy` is the shear
stress on interior cell corners. `extra` is the Schur term from the
semi-implicit thickness scheme (Arthern et al. 2015).
"""
@kernel function _op_h_stresses!(r_xx, r_yy, extra, r_xy, u, v, u_idx, v_idx, hu, hv,
                                 mask_u, mask_v, mask_c, D_h, D_imp, D_c,
                                 dx_inv, dy_inv, vecSampled)
    i, j = @index(Global, NTuple)
    @inbounds begin
        z = zero(eltype(u))
        u_ij = velocity_at(u, u_idx, i, j, vecSampled, z)
        u_ip = velocity_at(u, u_idx, i + 1, j, vecSampled, z)
        v_ij = velocity_at(v, v_idx, i, j, vecSampled, z)
        v_jp = velocity_at(v, v_idx, i, j + 1, vecSampled, z)
        dudx = (u_ip - u_ij) * dx_inv
        dvdy = (v_jp - v_ij) * dy_inv
        D = D_h[i, j]
        r_xx[i, j] = -2 * D * (2 * dudx + dvdy)
        r_yy[i, j] = -2 * D * (dudx + 2 * dvdy)

        qx_l = mask_u[i, j] ? hu[i, j] * u_ij : z
        qx_r = mask_u[i + 1, j] ? hu[i + 1, j] * u_ip : z
        qy_b = mask_v[i, j] ? hv[i, j] * v_ij : z
        qy_t = mask_v[i, j + 1] ? hv[i, j + 1] * v_jp : z
        extra[i, j] = D_imp[i, j] * ((qx_r - qx_l) * dx_inv + (qy_t - qy_b) * dy_inv)

        nxc = size(r_xy, 1)
        nyc = size(r_xy, 2)
        if i <= nxc && j <= nyc
            u_ipjp = velocity_at(u, u_idx, i + 1, j + 1, vecSampled, z)
            v_ipjp = velocity_at(v, v_idx, i + 1, j + 1, vecSampled, z)
            dudy = (u_ipjp - u_ip) * dy_inv
            dvdx = (v_ipjp - v_jp) * dx_inv
            r_xy[i, j] = mask_c[i, j] ? -D_c[i, j] * (dudy + dvdx) : zero(eltype(r_xy))
        end
    end
end

"""
    _op_force_u!(out_u, r_xx, r_xy, extra, u, u_idx, hu, D_u, dx_inv, dy_inv, vecSampled)

Write the x-force at each inner u-point from neighbouring stresses, basal drag,
and the Schur term (Arthern et al. 2015).

Skip the point when `u_idx[i, j]` is `0`, which means ocean or a fixed
boundary. Otherwise write the force into `out_u[u_idx[i, j]]`, the same
slot as that point's velocity in the short inner list.
"""
@kernel function _op_force_u!(out_u, r_xx, r_xy, extra, u, u_idx, hu, D_u, dx_inv, dy_inv, vecSampled)
    i, j = @index(Global, NTuple)
    @inbounds begin
        k = u_idx[i, j]
        if k != 0
            nxh = size(r_xx, 1)
            nxc = size(r_xy, 1)
            nyc = size(r_xy, 2)
            z = zero(eltype(out_u))

            r_left = i > 1 ? r_xx[i - 1, j] : z
            r_right = i <= nxh ? r_xx[i, j] : z
            d_rxx_dx = (r_right - r_left) * dx_inv

            if i > 1 && i <= nxc + 1
                val_bot = j > 1 ? r_xy[i - 1, j - 1] : z
                val_top = j <= nyc ? r_xy[i - 1, j] : z
                d_rxy_dy = (val_top - val_bot) * dy_inv
            else
                d_rxy_dy = z
            end

            e_left = i > 1 ? extra[i - 1, j] : z
            e_right = i <= nxh ? extra[i, j] : z
            d_extra_dx = (e_right - e_left) * dx_inv

            u_ij = velocity_at(u, u_idx, i, j, vecSampled, z)
            taubx = -D_u[i, j] * u_ij
            out_u[k] = d_rxx_dx + d_rxy_dy - taubx - hu[i, j] * d_extra_dx
        end
    end
end

"""
    _op_force_v!(out_v, r_yy, r_xy, extra, v, v_idx, hv, D_v, dx_inv, dy_inv, vecSampled)

Write the y-force at each inner v-point from neighbouring stresses, basal drag,
and the Schur term (Arthern et al. 2015).

Skip the point when `v_idx[i, j]` is `0`, which means ocean or a fixed
boundary. Otherwise write the force into `out_v[v_idx[i, j]]`, the same
slot as that point's velocity in the short inner list.
"""
@kernel function _op_force_v!(out_v, r_yy, r_xy, extra, v, v_idx, hv, D_v, dx_inv, dy_inv, vecSampled)
    i, j = @index(Global, NTuple)
    @inbounds begin
        k = v_idx[i, j]
        if k != 0
            nyh = size(r_yy, 2)
            nxc = size(r_xy, 1)
            nyc = size(r_xy, 2)
            z = zero(eltype(out_v))

            r_bot = j > 1 ? r_yy[i, j - 1] : z
            r_top = j <= nyh ? r_yy[i, j] : z
            d_ryy_dy = (r_top - r_bot) * dy_inv

            if j > 1 && j <= nyc + 1
                val_left = i > 1 ? r_xy[i - 1, j - 1] : z
                val_right = i <= nxc ? r_xy[i, j - 1] : z
                d_rxy_dx = (val_right - val_left) * dx_inv
            else
                d_rxy_dx = z
            end

            e_bot = j > 1 ? extra[i, j - 1] : z
            e_top = j <= nyh ? extra[i, j] : z
            d_extra_dy = (e_top - e_bot) * dy_inv

            v_ij = velocity_at(v, v_idx, i, j, vecSampled, z)
            tauby = -D_v[i, j] * v_ij
            out_v[k] = d_ryy_dy + d_rxy_dx - tauby - hv[i, j] * d_extra_dy
        end
    end
end

"""
    _op_diag_u!(op_diag_u, inner_indices, D_h, D_c, D_u, D_imp, hu, mask_u, mask_c, dx_inv, dy_inv)

The Gauss-Seidel weight at each inner u-point is the force on this point
if its velocity is 1 and every other velocity is 0. That is the same
physics as `_op_h_stresses!` plus `_op_force_u!`, but only the diagonal
entry.

`inner_indices[k]` is the 2D grid location of the k-th inner point.
"""
@kernel function _op_diag_u!(op_diag_u, inner_indices, D_h, D_c, D_u, D_imp,
                             hu, mask_u, mask_c, dx_inv, dy_inv)
    k = @index(Global, Linear)
    @inbounds begin
        nxu = size(D_u, 1)
        nxh = size(D_h, 1)
        nxc = size(D_c, 1)
        nyc = size(D_c, 2)
        idx = inner_indices[k]
        i = (idx - 1) % nxu + 1
        j = (idx - 1) ÷ nxu + 1
        z = zero(eltype(op_diag_u))
        dx2 = dx_inv * dx_inv
        dy2 = dy_inv * dy_inv

        d_h = z
        if i <= nxh
            d_h += D_h[i, j]
        end
        if i > 1
            d_h += D_h[i - 1, j]
        end
        diag = (4 * dx2) * d_h

        if i > 1 && i <= nxc + 1
            if j <= nyc && mask_c[i - 1, j]
                diag += D_c[i - 1, j] * dy2
            end
            if j > 1 && mask_c[i - 1, j - 1]
                diag += D_c[i - 1, j - 1] * dy2
            end
        end

        diag += D_u[i, j]

        if mask_u[i, j]
            d_imp = z
            if i <= nxh
                d_imp += D_imp[i, j]
            end
            if i > 1
                d_imp += D_imp[i - 1, j]
            end
            h = hu[i, j]
            diag += (h * h * dx2) * d_imp
        end

        op_diag_u[k] = diag
    end
end

"""
    _op_diag_v!(op_diag_v, inner_indices, D_h, D_c, D_v, D_imp, hv, mask_v, mask_c, dx_inv, dy_inv)

The Gauss-Seidel weight at each inner v-point is the force on this point
if its velocity is 1 and every other velocity is 0. That is the same
physics as `_op_h_stresses!` plus `_op_force_v!`, but only the diagonal
entry.

`inner_indices[k]` is the 2D grid location of the k-th inner point.
"""
@kernel function _op_diag_v!(op_diag_v, inner_indices, D_h, D_c, D_v, D_imp,
                             hv, mask_v, mask_c, dx_inv, dy_inv)
    k = @index(Global, Linear)
    @inbounds begin
        nxv = size(D_v, 1)
        nyh = size(D_h, 2)
        nxc = size(D_c, 1)
        nyc = size(D_c, 2)
        idx = inner_indices[k]
        i = (idx - 1) % nxv + 1
        j = (idx - 1) ÷ nxv + 1
        z = zero(eltype(op_diag_v))
        dx2 = dx_inv * dx_inv
        dy2 = dy_inv * dy_inv

        d_h = z
        if j <= nyh
            d_h += D_h[i, j]
        end
        if j > 1
            d_h += D_h[i, j - 1]
        end
        diag = (4 * dy2) * d_h

        if j > 1 && j <= nyc + 1
            if i <= nxc && mask_c[i, j - 1]
                diag += D_c[i, j - 1] * dx2
            end
            if i > 1 && mask_c[i - 1, j - 1]
                diag += D_c[i - 1, j - 1] * dx2
            end
        end

        diag += D_v[i, j]

        if mask_v[i, j]
            d_imp = z
            if j <= nyh
                d_imp += D_imp[i, j]
            end
            if j > 1
                d_imp += D_imp[i, j - 1]
            end
            h = hv[i, j]
            diag += (h * h * dy2) * d_imp
        end

        op_diag_v[k] = diag
    end
end

end # module
