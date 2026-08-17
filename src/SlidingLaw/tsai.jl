export TsaiSlidingLaw

struct TsaiSlidingLaw{T <: Real, C <: Union{T,Array{T,2}}, W <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    coulomb_coefficient :: C
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
end

"""
TsaiSlidingLaw(; <kwargs>)


Keyword arguments
=================
- coulomb_coefficient  : Coulomb friction coefficient
- drag_coefficient     : Weertman friction coefficients [Pa (yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
"""

function TsaiSlidingLaw(;
                        coulomb_coefficient = 0.5,
                        drag_coefficient = 1.0e4,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5)

    return TsaiSlidingLaw(
                            coulomb_coefficient,
                            drag_coefficient,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::TsaiSlidingLaw, model::AbstractModel)

use Tsai sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::TsaiSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    update_drag_coefficient!(model)
    gh.β .= (min.((sliding_law.coulomb_coefficient .* gh.effective_pressure) ./ ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ),
                gh.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)))
    return model
end

function reconstruct_on_grid(sliding_law::TsaiSlidingLaw, grid::Grid)
    return TsaiSlidingLaw(
        field_on_grid(sliding_law.coulomb_coefficient, grid; name = "Coulomb Coefficient"),
        field_on_grid(sliding_law.drag_coefficient, grid; name = "Drag Coefficient"),
        sliding_law.weertman_m,
        sliding_law.reg_speed)
end

function reconstruct_on_subdomain(sliding_law::TsaiSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    return TsaiSlidingLaw(
          spatial_on_subdomain(sliding_law.coulomb_coefficient, grid, subdomain),
          spatial_on_subdomain(sliding_law.drag_coefficient, grid, subdomain),
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end