export BuddSlidingLaw

struct BuddSlidingLaw{T <: Real, W <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
    budd_q :: T
end

"""
BuddSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_coefficient     : Budd friction coefficients [(yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
- budd_q               : Budd exponent
"""

function BuddSlidingLaw(;
                        drag_coefficient = 0.117, # 0.117 (yr/m)^(1/3) = 37.01 (s/m)^(1/3)
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5,
                        budd_q = 1.0)

    return BuddSlidingLaw(
                            drag_coefficient,
                            weertman_m,
                            reg_speed,
                            budd_q)
end

"""
            update_β_using_sliding_law!(sliding_law::BuddSlidingLaw, model::AbstractModel)

use Budd sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::BuddSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    gh.β .= sliding_law.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0) .* gh.effective_pressure.^(sliding_law.budd_q)
    return model
end

function reconstruct_on_grid(sliding_law::BuddSlidingLaw, grid::Grid)
    return BuddSlidingLaw(
        field_on_grid(sliding_law.drag_coefficient, grid; name = "Drag Coefficient"),
        sliding_law.weertman_m,
        sliding_law.reg_speed,
        sliding_law.budd_q)
end

function reconstruct_on_subdomain(sliding_law::BuddSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    return BuddSlidingLaw(
          spatial_on_subdomain(sliding_law.drag_coefficient, grid, subdomain),
          sliding_law.weertman_m,
          sliding_law.reg_speed,
          sliding_law.budd_q)
end