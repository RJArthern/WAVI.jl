export WeertmanSlidingLaw

struct WeertmanSlidingLaw{T <: Real, W <: Union{T,AbstractArray{T,2}}} <: AbstractSlidingLaw
    drag_coefficient :: W
    weertman_m :: T
    reg_speed :: T
end

"""
WeertmanSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_coefficient     : Weertman friction coefficients [Pa (yr/m)^(1/weertman_m)]
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
"""

function WeertmanSlidingLaw(; 
                        drag_coefficient = 1.0e4,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5) 
                        
    return WeertmanSlidingLaw(
                            drag_coefficient,
                            weertman_m,
                            reg_speed)
end

"""
            update_β_using_sliding_law!(sliding_law::WeertmanSlidingLaw, model::AbstractModel)

use Weertman sliding law to calculate basal drag coefficient
"""

function update_β_using_sliding_law!(sliding_law::WeertmanSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    update_drag_coefficient!(model)
    gh.β .= gh.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)
    return model
end

"""
    reconstruct_on_grid(sliding_law::WeertmanSlidingLaw, grid)

Return a copy whose drag coefficient is a 2D array on `grid`.

A number (the usual constructor default) is expanded to every cell. An array
must already match the grid. Called when the model is built, including after
an MPI subdomain or thread tile has been cut.
"""
function reconstruct_on_grid(sliding_law::WeertmanSlidingLaw, grid::Grid)
    return WeertmanSlidingLaw(
        field_on_grid(sliding_law.drag_coefficient, grid; name = "Drag Coefficient"),
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end

"""
    reconstruct_on_subdomain(sliding_law::WeertmanSlidingLaw, grid, subdomain)

Return a copy whose drag coefficient is cut to this tile, if it is already a
full-grid array. A number is left as a number so `reconstruct_on_grid` can
expand it on the local grid rather than the full domain.
"""
function reconstruct_on_subdomain(sliding_law::WeertmanSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    return WeertmanSlidingLaw(
          spatial_on_subdomain(sliding_law.drag_coefficient, grid, subdomain),
          sliding_law.weertman_m,
          sliding_law.reg_speed)
end

"""
    on_architecture(arch, sliding_law::WeertmanSlidingLaw)

Copy the drag coefficient onto `arch` when it is a dense array.
"""
function on_architecture(arch::AbstractArchitecture, sliding_law::WeertmanSlidingLaw)
    return WeertmanSlidingLaw(
        on_architecture(arch, sliding_law.drag_coefficient),
        sliding_law.weertman_m,
        sliding_law.reg_speed,
    )
end