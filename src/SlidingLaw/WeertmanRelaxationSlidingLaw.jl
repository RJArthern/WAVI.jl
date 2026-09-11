export WeertmanRelaxationSlidingLaw 

struct WeertmanRelaxationSlidingLaw{T <: Real, W <: Union{T,Array{T,2}},H <: Union{T,Array{T,2}}} <: AbstractSlidingLaw
    drag_coefficient :: W
    reference_thickness :: H
    weertman_m :: T
    reg_speed :: T
    convergence_timescale :: T
end

"""
WeertmanRelaxationSlidingLaw(; <kwargs>)


Keyword arguments
=================
- drag_coefficient     : Weertman friction coefficients [Pa (yr/m)^(1/weertman_m)]
- reference_thickness  : Target thickness field for relaxation (m) 
- weertman_m           : Weertman exponent
- reg_speed            : regularization speed, used to prevent bed speed going to zero
- convergence timescale : parameter to control rate of adjustment of drag coefficient (yr)
"""

function WeertmanRelaxationSlidingLaw(; 
                        drag_coefficient = 1.0e4,
                        reference_thickness = 2000.0,
                        weertman_m  = 3.0,
                        reg_speed = 1.0e-5,
                        convergence_timescale = 100.0) 
                        
    return WeertmanRelaxationSlidingLaw(
                            drag_coefficient,
                            reference_thickness,
                            weertman_m,
                            reg_speed,
                            convergence_timescale)
end

"""
            update_β_using_sliding_law!(sliding_law::WeertmanRelaxationSlidingLaw, model::AbstractModel)

use Weertman relaxation sliding law to calculate basal drag coefficient and update Weertman C coefficient to drive model towards a reference thickness.
"""

function update_β_using_sliding_law!(sliding_law::WeertmanRelaxationSlidingLaw, model::AbstractModel)
    @unpack gh=model.fields
    @unpack params=model
    update_drag_coefficient!(model)
    gh.β .= gh.drag_coefficient .* ( sqrt.(gh.bed_speed.^2 .+  sliding_law.reg_speed^2 ) ).^(1.0/sliding_law.weertman_m - 1.0)
    
    convergence_exponent = params.dt./sliding_law.convergence_timescale
    
    change_mask = gh.grounded_fraction .> 0
    sliding_law.drag_coefficient[change_mask] .= sliding_law.drag_coefficient[change_mask].*
           (sliding_law.reference_thickness[change_mask]./gh.h[change_mask]).^convergence_exponent

    return model
end

function reconstruct_on_grid(sliding_law::WeertmanRelaxationSlidingLaw, grid::Grid)
    return WeertmanRelaxationSlidingLaw(
        isa(sliding_law.drag_coefficient,Number) ? sliding_law.drag_coefficient*ones(grid.nx,grid.ny) : 
        size(sliding_law.drag_coefficient) == (grid.nx,grid.ny) ? sliding_law.drag_coefficient :
        throw(DimensionMismatch("Drag Coefficient does not match grid size")),
        isa(sliding_law.reference_thickness,Number) ? sliding_law.reference_thickness*ones(grid.nx,grid.ny) : 
        size(sliding_law.reference_thickness) == (grid.nx,grid.ny) ? sliding_law.reference_thickness :
        throw(DimensionMismatch("Reference thickness does not match grid size")),
          sliding_law.weertman_m,
          sliding_law.reg_speed,
          sliding_law.convergence_timescale
          )
end

function reconstruct_on_subdomain(sliding_law::WeertmanRelaxationSlidingLaw, grid::Grid, subdomain::NTuple{4,<: Integer})
    
    x_start,x_end,y_start,y_end = subdomain

    return WeertmanRelaxationSlidingLaw(
          sliding_law.drag_coefficient[x_start:x_end, y_start:y_end],
          sliding_law.reference_thickness[x_start:x_end, y_start:y_end],
          sliding_law.weertman_m,
          sliding_law.reg_speed,
          sliding_law.convergence_timescale)
end