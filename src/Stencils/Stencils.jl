module Stencils

export _diff_x!,
    _diff_y!,
    _diff_xT!,
    _diff_yT!,
    _avg_x!,
    _avg_y!,
    _apply_mask!,
    _scale!,
    _gather!,
    _scatter!,
    launch!

using KernelAbstractions: KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

"""
    launch!(kernel!, args...; ndrange, backend = KA.CPU())

Wrapper to launch a KernelAbstractions `kernel!` on the specified `backend`.
Defaults to `KA.CPU()` if no backend is provided. Ensures synchronisation
after the kernel is launched.
"""
function launch!(kernel!, args...; ndrange, backend = KA.CPU())
    kernel!(backend)(args...; ndrange = ndrange)
    KA.synchronize(backend)
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
"""
@kernel function _diff_xT!(out, inp, dx_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i - 1, j] - inp[i, j]) * dx_inv
end

"""
    _diff_yT!(out, inp, dy_inv)

Compute backward (transpose) finite difference of `inp` in y-direction.
Equivalent to `out = ∂yᵀ * inp` in sparse matrix notation.
"""
@kernel function _diff_yT!(out, inp, dy_inv)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = (inp[i, j - 1] - inp[i, j]) * dy_inv
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
