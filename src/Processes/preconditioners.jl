using LinearAlgebra
using Parameters

import WAVI: AbstractModel, AbstractSpec
using WAVI.Wavelets
using WAVI.Utilities

export update_preconditioner!, precondition!

function update_preconditioner!(model::AbstractModel{T,N,<:AbstractSpec}) where {T,N}
    return model
end

function precondition!(model::AbstractModel{T,N,<:AbstractSpec}) where {T,N}
    @unpack solver_params=model
    s = stencil_scratch!(model)

    x = get_start_guess(model)
    op = get_op(model)
    b = get_rhs(model)
    resid = s.picard_resid
    get_resid!(resid, x, op, b)
    set_residual!(model, resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    correction = s.picard_correction
    fill!(correction, zero(eltype(correction)))

    if ! converged
      p = get_preconditioner(model, op)
      apply_preconditioning!(correction, p, resid)
      correction_coarse = get_correction_coarse(p)
      set_correction_coarse!(model,correction_coarse)
    end
    x .= x .+ correction
    set_velocities!(model,x)

    return converged, rel_resid
end
