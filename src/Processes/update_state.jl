export update_state!, update_velocities_on_h_grid!, update_state_novelocity!

using Parameters

using WAVI: AbstractModel
using WAVI.KroneckerProducts
using WAVI.MeltRates
using WAVI.SurfaceMassBalance
using WAVI.Fracture
using WAVI.SlidingLaw
using WAVI.BasalHydrology
using WAVI.ThermoDynamics
using WAVI.Time
using WAVI.Utilities
using WAVI.Wavelets
using LinearAlgebra: Diagonal

# Sparse samp/spread and Kronecker maps (cent, ∂x, ∂y) stay on the host, and
# CuArray does not support boolean indexing such as dest[mask] .= packed.
# Copy dense fields to Array, apply those operators, then write back.
function _copy_masked!(dest, mask, packed)
    dest_h = _host(dest)
    dest_h[_host(mask)] .= packed
    dest_h === dest || copyto!(dest, dest_h)
    return dest
end

"""
update_state!(model::AbstractModel, clock)

Update the model to the current time dependent situation
"""
function update_state!(model::AbstractModel, clock::Clock) 
    @debug "Updating state at $(clock)"
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model, clock)
    update_shelf_basal_melt!(model, clock)
    update_thermodynamics_basal_melt!(model)
    update_glen_b!(model)
    update_dsdh!(model)
    update_basal_hydrology!(model; update_basal_water_thickness=true)
    update_model_velocities!(model)
    update_velocities_on_h_grid!(model)
    update_surf_speed!(model)
    update_strain_history!(model)
    update_dhdt!(model)
    update_model_wavelets!(model)
    update_surface_velocities_on_uv_grid!(model)
    return nothing
end

"""
update_state!(model::AbstractModel)

Update the model to the current time-indepdent situation
"""
function update_state!(model::AbstractModel) 
    @debug "Resolving state without time"
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model, Clock())
    update_shelf_basal_melt!(model, Clock())
    update_thermodynamics_basal_melt!(model)
    update_glen_b!(model)
    update_dsdh!(model)
    update_basal_hydrology!(model; update_basal_water_thickness=true)
    update_model_velocities!(model)
    update_velocities_on_h_grid!(model)
    update_surf_speed!(model)
    update_strain_history!(model)
    update_dhdt!(model)
    update_model_wavelets!(model)
    update_surface_velocities_on_uv_grid!(model)
    return nothing
end

"""
update_state_novelocity!(model::AbstractModel, clock)

Update the model to the current time dependent situation, without updating the ice velocities
"""
function update_state_novelocity!(model, clock)
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model, clock)
    update_shelf_basal_melt!(model, clock)
    update_thermodynamics_basal_melt!(model)
    update_glen_b!(model)
    update_dsdh!(model)
    update_basal_hydrology!(model; update_basal_water_thickness=true)
    update_velocities_on_h_grid!(model)
    update_surf_speed!(model)
    update_strain_history!(model)
    update_dhdt!(model)
    update_model_wavelets!(model)
    update_surface_velocities_on_uv_grid!(model)
    return nothing
end



"""
    update_surface_elevation!(model::AbstractModel)

Adjust surface elevation to hydrostatic equilibrium.
"""
function update_surface_elevation!(model::AbstractModel)
    @unpack params=model
    @unpack gh=model.fields
    sea = params.sea_level_wrt_geoid
    ρ = params.density_ice / params.density_ocean
    @. gh.s = ifelse(gh.mask, max(gh.b + gh.h, sea + gh.h * (1 - ρ)), gh.s)
    return model
end

"""
    update_geometry_on_uv_grids!(model::AbstractModel)

Interpolate thickness and surface elvation from h-grid to u- and v-grids.

"""
function update_geometry_on_uv_grids!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    # samp and centᵀ are host sparse/Kronecker operators; see _host / _copy_masked! above.
    h = _host(gh.h)
    s = _host(gh.s)
    crop = Diagonal(_host(gh.crop.diag))
    onesvec = ones(T, length(h))
    denu = gu.samp * (gu.centᵀ * (crop * onesvec))
    denv = gv.samp * (gv.centᵀ * (crop * onesvec))
    _copy_masked!(gu.h, gu.mask, (gu.samp * (gu.centᵀ * (crop * vec(h)))) ./ denu)
    _copy_masked!(gu.s, gu.mask, (gu.samp * (gu.centᵀ * (crop * vec(s)))) ./ denu)
    _copy_masked!(gv.h, gv.mask, (gv.samp * (gv.centᵀ * (crop * vec(h)))) ./ denv)
    _copy_masked!(gv.s, gv.mask, (gv.samp * (gv.centᵀ * (crop * vec(s)))) ./ denv)
    return model
end

"""
    update_height_above_floatation!(model::AbstractModel)

Update height above floatation. Zero value is used to define location of grounding line.
"""
function update_height_above_floatation!(model::AbstractModel)
    @unpack params=model
    @unpack gh=model.fields
    # Do not broadcast `params`: it holds host Matrix fields and cannot enter a GPU kernel.
    ρ = params.density_ocean / params.density_ice
    sea = params.sea_level_wrt_geoid
    @. gh.haf = gh.h - ρ * (sea - gh.b)
    return model
end

"""
    update_grounded_fraction_on_huv_grids!(model::AbstractModel)

Update grounded area fraction on h-, u-, and v-grids for use in subgrid parameterisation.
"""
function update_grounded_fraction_on_huv_grids!(model::AbstractModel)
    @unpack gh,gu,gv = model.fields
    # pos_fraction uses host boolean indexing; copyto! writes back to device arrays.
    (gfh,gfu,gfv)=pos_fraction(_host(gh.haf);mask=_host(gh.mask))
    copyto!(gh.grounded_fraction, gfh)
    copyto!(gu.grounded_fraction, gfu)
    copyto!(gv.grounded_fraction, gfv)
    return model
end

"""
    update_shelf_basal_melt!(model::AbstractModel)

Update the basal melt rate under ice shelves.
"""
function update_shelf_basal_melt!(model::AbstractModel, clock)
    update_shelf_melt_rate!(model.shelf_melt_rate, model.fields, model.grid, clock)
    return model
end

"""
    update_thermodynamics!(model::AbstractModel)

Update the ice temperature and grounded melt rate according to the chosen thermodynamics model.
The specific function lives in the corresponding thermodynamics file.
"""
function update_thermodynamics_basal_melt!(model::AbstractModel)
    update_ice_temperature_and_basal_melt_rate!(model.thermo_dynamics,model)
    return model
end

"""
    update_glen_b!(model::AbstractModel)

Update stiffness parameter B in Glen flow law.
"""
function update_glen_b!(model::AbstractModel)
    @unpack g3d = model.fields
    @unpack params = model
    fill_glen_b!(
        g3d.glen_b,
        g3d.θ,
        g3d.Φ,
        params.glen_a_ref,
        params.glen_n,
        params.glen_a_activation_energy,
        params.glen_temperature_ref,
        params.gas_const,
    )
    return model
end



"""
    update_dsdh!(model::AbstractModel)

Compute change of surface elevation per unit thickness change, accounting for hydrostatic adjustment.
"""
function update_dsdh!(model::AbstractModel)
    @unpack gh,gu,gv=model.fields
    @unpack params = model
    gh.dsdh .= (1.0 - params.density_ice./params.density_ocean) .+
           (params.density_ice./params.density_ocean).*gh.grounded_fraction;
    return model
end

"""
    update_basal_hydrology!(model::AbstractModel; update_basal_water_thickness::Bool = true)

Update the basal water thickness and effective pressure according to the chosen basal hydrology model.
The specific function lives in the corresponding basal hydrology file.
"""
function update_basal_hydrology!(model::AbstractModel; update_basal_water_thickness::Bool = true)
    update_basal_water_thickness_effective_pressure!(model.basal_hydrology,model; update_basal_water_thickness=update_basal_water_thickness)
    return model
end

"""
    update_model_velocity!(model::AbstractModel)

Wrapper function for that which updates the model velocities on the u, v grids (update_velocities in separate file)
"""
function update_model_velocities!(model::AbstractModel)
    update_velocities!(model)
    return model
end

"""
    update_velocities_on_h_grid!(model::AbstractModel)

Update the velocities (depth averaged, surface and bed) on the h grid 
"""
function update_velocities_on_h_grid!(model::AbstractModel{T,N,S}) where {T,N,S<:AbstractSpec}
    @unpack gh,gu,gv = model.fields
    #depth averaged velocities
    copyto!(gh.u, reshape(gu.cent * vec(_host(gu.u)), size(gh.u)))
    copyto!(gh.v, reshape(gv.cent * vec(_host(gv.v)), size(gh.v)))

    #bed velocities
    gh.ub .= gh.u ./ (1 .+ (gh.β .* gh.quad_f2))
    gh.vb .= gh.v ./ (1 .+ (gh.β .* gh.quad_f2))

    #surface velocities
    gh.us .= gh.ub .* (1 .+ (gh.β .* gh.quad_f1))
    gh.vs .= gh.vb .* (1 .+ (gh.β .* gh.quad_f1))
    return model
end

"""
    update_surf_speed!(model::AbstractModel)

Find the sliding speed on the h-grid using the speed components.
"""
function update_surf_speed!(model::AbstractModel)
    @unpack gh=model.fields
    gh.surf_speed  .= sqrt.(gh.us.^2 .+gh.vs.^2);
    return model
end

"""
    update_dhdt!(model::AbstractModel)

Evaluate rate of change of thickness using mass conservation.
"""
function update_dhdt!(model::AbstractModel)
    @unpack gh,gu,gv=model.fields
    crop_u = Diagonal(_host(gu.crop.diag))
    crop_v = Diagonal(_host(gv.crop.diag))
    packed = gh.samp * (
        vec(_host(gh.accumulation)) .- vec(_host(gh.basal_melt)) .- (
            (gu.∂x * (crop_u * (vec(_host(gu.h)) .* vec(_host(gu.u))))) .+
            (gv.∂y * (crop_v * (vec(_host(gv.h)) .* vec(_host(gv.v)))))
        )
    )
    _copy_masked!(gh.dhdt, gh.mask, packed)
    return model
end

""" 
    update_model_wavelets(model::AbstractModel)

Wrapper function for that which updates the model wavelets
"""
function update_model_wavelets!(model::AbstractModel)
    update_wavelets!(model)
    return model
end

function update_surface_velocities_on_uv_grid!(model)
    @unpack gh,gu,gv = model.fields
    #surface  velocities
    crop_h = Diagonal(_host(gh.crop.diag))
    crop_u = Diagonal(_host(gu.crop.diag))
    crop_v = Diagonal(_host(gv.crop.diag))
    copyto!(gu.us, reshape(crop_u * (gu.centᵀ * (crop_h * vec(_host(gh.us)))), size(gu.us)))
    copyto!(gv.vs, reshape(crop_v * (gv.centᵀ * (crop_h * vec(_host(gh.vs)))), size(gv.vs)))
    return model
end

