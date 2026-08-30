export CoulombSlidingLaw

struct CoulombSlidingLaw{T <: Real, C <: Union{T,AbstractArray{T,2}}} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    reg_speed :: T
end

"""
CoulombSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient    : Coulomb friction coefficient
- reg_speed : regularization speed, used to prevent bed speed going to zero
"""

function CoulombSlidingLaw(;
                        coulomb_coefficient = 0.5,
                        reg_speed=1.0e-5)

    return CoulombSlidingLaw(
                            coulomb_coefficient,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel)

use Coulomb sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::CoulombSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= (sliding_law.coulomb_coefficient .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) )
    return model
end


function reconstruct_on_grid(sliding_law::CoulombSlidingLaw, grid::Grid)
    return CoulombSlidingLaw(
        field_on_grid(sliding_law.coulomb_coefficient, grid; name = "Coulomb Coefficient"),
        sliding_law.reg_speed)
end

function reconstruct_on_subdomain(sliding_law::CoulombSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    return CoulombSlidingLaw(
          spatial_on_subdomain(sliding_law.coulomb_coefficient, grid, subdomain),
          sliding_law.reg_speed)
end

function on_architecture(arch::AbstractArchitecture, sliding_law::CoulombSlidingLaw)
    return CoulombSlidingLaw(
        on_architecture(arch, sliding_law.coulomb_coefficient),
        sliding_law.reg_speed,
    )
end