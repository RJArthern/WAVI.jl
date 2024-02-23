
"""
update_state!(model::AbstractModel, clock)

Update the model to the current time dependent situation
"""
function update_state!(model, clock)
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model)
    update_basal_melt!(model, clock)
    update_weertman_c!(model)
    update_dsdh!(model)
    update_model_velocities!(model)
    update_velocities_on_h_grid!(model)
    update_dhdt!(model)
    update_model_wavelets!(model)
    return nothing
end

"""
update_state!(model::AbstractModel)

Update the model to the current time-indepdent situation
"""
function update_state!(model)
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model)
    update_basal_melt!(model, WAVI.Clock())
    update_weertman_c!(model)
    update_dsdh!(model)
    update_model_velocities!(model)
    update_velocities_on_h_grid!(model)
    update_dhdt!(model)
    update_model_wavelets!(model)
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
    
    if params.step_haf == false
    if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny
 
    haf_nochild=gh.haf[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]  

    end
    end
    
    gh.haf .= height_above_floatation.(gh.h,gh.b,Ref(params))
    
    if params.step_haf == false
    if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny
    gh.haf[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]=haf_nochild[:,:]
    end
    end
    
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
    update_accumulation_rate!(model::AbstractModel)

Update the accumulation rate.
"""
function update_accumulation_rate!(model::AbstractModel)
    @unpack params = model
    @unpack gh=model.fields
    gh.accumulation .= params.accumulation_rate
    return model
end


"""
    update_basal_melt!(model::AbstractModel)

Update the basal melt rate.
"""
function update_basal_melt!(model::AbstractModel, clock)
    update_melt_rate!(model.melt_rate, model.fields, model.grid, clock)
    return model
end

"""
    update_weertman_c!(model::AbstractModel)

Update coefficient used in the sliding law to account for migration of grounding line.
"""
function update_weertman_c!(model::AbstractModel)
    @unpack gh=model.fields
    @unpack params=model
    gh.weertman_c .= params.weertman_c .* gh.grounded_fraction
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
function update_velocities_on_h_grid!(model)
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
