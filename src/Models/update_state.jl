
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
    update_weertman_c!(model,clock)
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
    update_weertman_c!(model, WAVI.Clock())
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
function update_geometry_on_uv_grids!(model::AbstractModel)
    @unpack gh,gu,gv,gc=model.fields
    onesvec=ones(gh.nxh*gh.nyh)
    gu.h[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.h[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    gu.s[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.s[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    gv.h[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.h[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
    gv.s[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.s[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
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
    @unpack basal_melt,grounded_fraction=model.fields.gh
    @unpack params,grid=model
    update_melt_rate!(model.melt_rate, model.fields, model.grid, clock)
    #modify for tidal signal
    if params.tidal_melting
    (nx,ny) = size(model.grid.xxh)
    #A = get_normalized_tidal_amplitude(clock.time, params.tidal_daily_timescale, params.tidal_hourly_timescale)
    #println(A)
    for ix = 1:nx
        for iy = 1:ny
            if (model.fields.gh.grounded_fraction[ix,iy] > 0)
                #get the tidal amplitude at the current time
                A = get_normalized_tidal_amplitude(clock.time, params.tidal_daily_timescale, params.tidal_hourly_timescale)
                L = A*params.tidal_lengthscale
                if (ix,iy) == (1,1)
                    show(A)
                    show("   ")
                end

                # find the distance to the closest fully floating point and the melt there
                dp = sqrt.((grid.xxh .- grid.xxh[ix,iy]).^2 .+ (grid.yyh .- grid.yyh[ix,iy]).^2);
                dp[grounded_fraction .> 0] .= 1e16; #set grounded ice and partially floatign to large, so we never pick this up
                dgl = findmin(dp)[1]
                mgl = basal_melt[findmin(dp)[2]]
                # set the basal melt rate
                basal_melt[ix,iy] = mgl * exp.(-(dgl).^2 / 2 / L.^2)
            end

        end
    end
    end
    return model
end

"""
    update_weertman_c!(model::AbstractModel)

Update coefficient used in the sliding law to account for migration of grounding line.
"""
function update_weertman_c!(model::AbstractModel,clock)
    @unpack gh=model.fields
    @unpack params,grid=model
    if ~(params.tidal_drag) #standard drag formulation 
        if params.partial_cell_drag 
            gh.weertman_c .= params.weertman_c .* gh.grounded_fraction
        else #no partial cell drag
            grf = deepcopy(gh.grounded_fraction)
            grf[grf .> 0] .= 1      #set anything partially floating to grounded 
            gh.weertman_c  .=  params.weertman_c .* grf
        end

    else # tidal melting formultation 
        #get the tidal amplitude at the current time
        A = get_normalized_tidal_amplitude(clock.time, params.tidal_daily_timescale, params.tidal_hourly_timescale)
        L = A*params.tidal_lengthscale
        dgl = get_grounding_line_distance(gh.grounded_fraction,grid.xxh,grid.yyh)
        tidal_modification = (1 .- exp.(-(dgl).^2 / 2 / L.^2));
        gh.weertman_c .= params.weertman_c .* tidal_modification
    end
    
    return model
end

"""
    get_normalized_tidal_amplitude(clock)

Return the normalized tidal amplitude associated with the current time
"""
function get_normalized_tidal_amplitude(t, tidal_daily_timescale, tidal_hourly_timescale)
    tt = t*365.25 * 24 #convert to hours
    amp = cos(2*pi*tt / (14*24))*cos(2*pi*tt / 6.25); 
    amp = cos(2*pi*tt / (tidal_daily_timescale*24))*cos(2*pi*tt / tidal_hourly_timescale); 
    return 1/2*(1 + amp) #between 0 and 1
end

"""
    get_grounding_line_distance(grounded_fraction,xx,yy)

Return an array of distance to the nearest fully floating point (i.e. zero for the shelf)
"""

function  get_grounding_line_distance(grounded_fraction,xx,yy)
    (nx,ny) = size(grounded_fraction)
    idx = (grounded_fraction .== 0)
    M = zeros(nx,ny)
    for ix = 1:nx
        for iy = 1:ny
                dp = sqrt.((xx .- xx[ix,iy]).^2 .+ (yy .- yy[ix,iy]).^2);
                M[ix,iy] = minimum(dp[idx])

        end
    end
    return M

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