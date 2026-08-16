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
    _gather!,
    _scatter!,
    launch!

using KernelAbstractions: KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

"""
    launch!(kernel!, args...; ndrange, sync = true)

Wrapper to launch a KernelAbstractions `kernel!`.
The backend is taken from the first kernel argument via `KA.get_backend`.
Ensures synchronisation after the kernel is launched unless `sync = false`.
"""
function launch!(kernel!, args...; ndrange, sync = true)
    backend = KA.get_backend(first(args))
    kernel!(backend)(args...; ndrange = ndrange)
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
    _scale!(out, inp, diag)

Multiply each element of `inp` by the corresponding element of `diag`.
Equivalent to `out = Diagonal(diag) * inp`.

Multiplies every cell by a physical property, like viscosity.
"""
@kernel function _scale!(out, inp, diag)
    i = @index(Global, Linear)
    @inbounds out[i] = diag[i] * inp[i]
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

end # module
