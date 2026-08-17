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
    gh.s[gh.mask] .= max.(gh.b[gh.mask]+gh.h[gh.mask],
                          params.sea_level_wrt_geoid .+ gh.h[gh.mask]*(1-params.density_ice./params.density_ocean))
    return model
end

"""
    update_geometry_on_uv_grids!(model::AbstractModel)

Interpolate thickness and surface elvation from h-grid to u- and v-grids.

"""
function update_geometry_on_uv_grids!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    onesvec=ones(T,gh.nxh*gh.nyh)
    gu.h[gu.mask].=(gu.samp*(gu.centᵀ*(gh.crop*gh.h[:])))./(gu.samp*(gu.centᵀ*(gh.crop*onesvec)))
    gu.s[gu.mask].=(gu.samp*(gu.centᵀ*(gh.crop*gh.s[:])))./(gu.samp*(gu.centᵀ*(gh.crop*onesvec)))
    gv.h[gv.mask].=(gv.samp*(gv.centᵀ*(gh.crop*gh.h[:])))./(gv.samp*(gv.centᵀ*(gh.crop*onesvec)))
    gv.s[gv.mask].=(gv.samp*(gv.centᵀ*(gh.crop*gh.s[:])))./(gv.samp*(gv.centᵀ*(gh.crop*onesvec)))
    return model
end

"""
    update_height_above_floatation!(model::AbstractModel)

Update height above floatation. Zero value is used to define location of grounding line.
"""
function update_height_above_floatation!(model::AbstractModel)
    @unpack params=model
    @unpack gh=model.fields
    gh.haf .= height_above_floatation.(gh.h,gh.b,Ref(params))
    return model
end

"""
    update_grounded_fraction_on_huv_grids!(model::AbstractModel)

Update grounded area fraction on h-, u-, and v-grids for use in subgrid parameterisation.
"""
function update_grounded_fraction_on_huv_grids!(model::AbstractModel)
    @unpack gh,gu,gv = model.fields
    (gfh,gfu,gfv)=pos_fraction(gh.haf;mask=gh.mask)
    gh.grounded_fraction[:] .= gfh[:]
    gu.grounded_fraction[:] .= gfu[:]
    gv.grounded_fraction[:] .= gfv[:]
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
    gh.u[:] .= gu.cent*gu.u[:] #(gu.u[1:end-1,:] + gu.u[2:end,:])./2
    gh.v[:] .= gv.cent*gv.v[:] #(gv.v[:,1:end-1] + gv.v[:, 2:end])./2

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
    gh.dhdt[gh.mask].=gh.samp*(gh.accumulation[:] .- gh.basal_melt[:] .-
             (  (gu.∂x*(gu.crop*(gu.h[:].*gu.u[:]))) .+ (gv.∂y*(gv.crop*(gv.h[:].*gv.v[:]))) ) )
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
    gu.us[:].=gu.crop*(gu.centᵀ*gh.crop*(gh.us[:]))
    gv.vs[:].=gv.crop*(gv.centᵀ*gh.crop*(gh.vs[:]))
    return model
end

