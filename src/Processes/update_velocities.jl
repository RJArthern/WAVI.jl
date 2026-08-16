export update_velocities!, get_start_guess, get_op, get_rhs, get_resid, set_residual!

using LinearMaps

using WAVI
using WAVI.Parameters
using WAVI.Utilities
using WAVI.Stencils
using KernelAbstractions: KernelAbstractions as KA
using KernelAbstractions: @kernel, @index

"""
update_velocities!(model::AbstractModel)

Solve momentum equation to update the velocities, plus Picard iteration for non-linear rheology.

"""
function update_velocities!(model::AbstractModel{T,N}) where {T,N}
    @unpack solver_params = model

    update_preconditioner!(model)

    converged::Bool = false
    i_picard::Int64 = 0
    rel_resid = Inf
    while !converged && (i_picard < solver_params.maxiter_picard)

        i_picard = i_picard + 1
        @debug "Updating velocities, picard iteration $(i_picard)"
        inner_update!(model)
        
        converged, rel_resid = precondition!(model)

    end
    @debug "Solved momentum equation with residual $(round(rel_resid,sigdigits=3)) at iteration $(i_picard)"

    return model
end


function inner_update!(model::AbstractModel)
    update_shelf_strain_rate!(model)
    update_av_speed!(model)
    update_bed_speed!(model)
    update_basal_hydrology!(model; update_basal_water_thickness=false)
    update_β!(model)
    update_basal_drag!(model)
    update_glen_b!(model)
    update_damage!(model;store_strain_history = false)
    inner_update_viscosity!(model)
    update_av_viscosity!(model)
    update_quadrature_falpha!(model)
    update_βeff!(model)
    update_βeff_on_uv_grids!(model)
    update_rheological_operators!(model)
  #  update_surface_velocities_on_uv_grid!(model)
    return model
end

"""
    get_start_guess(model::AbstractModel)

 Return starting guess used to begin iterative solution of velocities.

"""
function get_start_guess(model::AbstractModel)
    @unpack gu,gv=model.fields
    @assert eltype(gu.u)==eltype(gv.v)
    x=[gu.samp_inner*gu.u[:];gv.samp_inner*gv.v[:]]
    return x
end


"""
    get_rhs(model::AbstractModel)

 Return right hand side vector of momentum equations.

"""
function get_rhs(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack params, solver_params = model
    grid = model.grid
    dx_inv = one(T) / grid.dx
    dy_inv = one(T) / grid.dy
    backend = KA.get_backend(gh.h)

    gu_inner_indices = findall(vec(gu.mask_inner))
    gv_inner_indices = findall(vec(gv.mask_inner))

    surf_crop = similar(gh.h)
    ones_crop = similar(gh.h)
    tmpu = similar(gu.u)
    tmpv = similar(gv.v)
    tmpui = similar(gu.u, gu.ni)
    tmpvi = similar(gv.v, gv.ni)

    rhs = zeros(T,gu.ni+gv.ni)
    f1 = zeros(T,gu.ni+gv.ni)
    f2 = zeros(T,gu.ni+gv.ni)
    f3 = zeros(T,gu.ni+gv.ni)
    sui = zeros(T,gu.ni)
    hui = zeros(T,gu.ni)
    dui = zeros(T,gu.ni)
    svi = zeros(T,gv.ni)
    hvi = zeros(T,gv.ni)
    dvi = zeros(T,gv.ni)

    @. surf_crop = gh.s + solver_params.super_implicitness * params.dt * gh.dsdh * (gh.accumulation - gh.basal_melt)
    launch!(_apply_mask!, surf_crop, gh.mask; ndrange = size(surf_crop), sync = false)
    fill!(ones_crop, one(T))
    launch!(_apply_mask!, ones_crop, gh.mask; ndrange = size(ones_crop), sync = false)
    KA.synchronize(backend)

    launch!(_diff_xT!, tmpu, surf_crop, dx_inv; ndrange = size(tmpu), sync = false)
    launch!(_diff_yT!, tmpv, surf_crop, dy_inv; ndrange = size(tmpv), sync = false)
    KA.synchronize(backend)
    @. tmpu = -tmpu
    @. tmpv = -tmpv
    launch!(_gather!, tmpui, tmpu, gu_inner_indices; ndrange = length(tmpui), sync = false)
    launch!(_gather!, tmpvi, tmpv, gv_inner_indices; ndrange = length(tmpvi), sync = false)
    KA.synchronize(backend)
    @. tmpui = (params.density_ice*params.g*gu.h[gu.mask_inner]).* tmpui
    @. tmpvi = (params.density_ice*params.g*gv.h[gv.mask_inner]).*tmpvi

    f1[1:gu.ni] .= tmpui
    f1[(gu.ni+1):(gu.ni+gv.ni)] .= tmpvi

    sui .= gu.s[gu.mask_inner]
    hui .= gu.h[gu.mask_inner]
    dui .= icedraft.(sui,hui,params.sea_level_wrt_geoid)
    launch!(_diff_xT!, tmpu, ones_crop, dx_inv; ndrange = size(tmpu))
    @. tmpu = -tmpu
    launch!(_gather!, tmpui, tmpu, gu_inner_indices; ndrange = length(tmpui))
    @. tmpui = tmpui * params.g*(0.5*params.density_ice*hui^2
                            - 0.5*params.density_ocean*dui^2
                            - params.density_ice*hui*sui)

    svi .= gv.s[gv.mask_inner]
    hvi .= gv.h[gv.mask_inner]
    dvi .= icedraft.(svi,hvi,params.sea_level_wrt_geoid)
    launch!(_diff_yT!, tmpv, ones_crop, dy_inv; ndrange = size(tmpv))
    @. tmpv = -tmpv
    launch!(_gather!, tmpvi, tmpv, gv_inner_indices; ndrange = length(tmpvi))
    @. tmpvi = tmpvi * params.g*(0.5*params.density_ice*hvi^2
                            - 0.5*params.density_ocean*dvi^2
                            - params.density_ice*hvi*svi)

    f2[1:gu.ni] .= tmpui
    f2[(gu.ni+1):(gu.ni+gv.ni)] .= tmpvi

 #   f2=[
 #       (0.5*params.density_ice*params.g*gu.h[gu.mask_inner].^2
 #       .- 0.5*params.density_ocean*params.g*(icedraft.(gu.s[gu.mask_inner],gu.h[gu.mask_inner],params.sea_level_wrt_geoid)).^2
 #       .- params.density_ice*params.g*gu.h[gu.mask_inner].*gu.s[gu.mask_inner]).*gu.samp_inner*(-gu.∂xᵀ*(gh.crop*onesvec))
 #       ;
 #       (0.5*params.density_ice*params.g*gv.h[gv.mask_inner].^2
 #       .- 0.5*params.density_ocean*params.g*(icedraft.(gv.s[gv.mask_inner],gv.h[gv.mask_inner],params.sea_level_wrt_geoid)).^2
 #       .- params.density_ice*params.g*gv.h[gv.mask_inner].*gv.s[gv.mask_inner]).*gv.samp_inner*(-gv.∂yᵀ*(gh.crop*onesvec))
 #       ]
    
    get_rhs_dirichlet!(f3,model)

    rhs .= f1 .+ f2 .+ f3

    return rhs
end


"""
    set_velocities!(model::AbstractModel,x)

Set velocities to particular values. Input vector x represents stacked u and v components at valid grid points.
"""
function set_velocities!(model::AbstractModel,x)
    @unpack gh,gu,gv,gc=model.fields
    @views gu.u[gu.mask_inner] .= x[1:gu.ni]
    @views gv.v[gv.mask_inner] .= x[(gu.ni+1):(gu.ni+gv.ni)]
    return model
end

"""
    update_shelf_strain_rate!(model::AbstractModel)

Find the effective strain rate for 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
function update_shelf_strain_rate!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    grid = model.grid
    dx_inv = one(T) / grid.dx
    dy_inv = one(T) / grid.dy
    backend = KA.get_backend(gh.h)

    u_crop = similar(gu.u)
    v_crop = similar(gv.v)
    dudx = similar(gh.h)
    dvdy = similar(gh.h)
    dudy_c = similar(gh.h, T, gc.nxc, gc.nyc)
    dvdx_c = similar(dudy_c)
    shear_c = similar(dudy_c)
    shear_h = similar(gh.h)

    copyto!(u_crop, gu.u)
    copyto!(v_crop, gv.v)
    launch!(_apply_mask!, u_crop, gu.mask; ndrange = size(u_crop), sync = false)
    launch!(_apply_mask!, v_crop, gv.mask; ndrange = size(v_crop), sync = false)
    KA.synchronize(backend)

    launch!(_diff_x!, dudx, u_crop, dx_inv; ndrange = size(dudx), sync = false)
    launch!(_diff_y!, dvdy, v_crop, dy_inv; ndrange = size(dvdy), sync = false)
    launch!(_diff_y_staggered!, dudy_c, u_crop, dy_inv; ndrange = size(dudy_c), sync = false)
    launch!(_diff_x_staggered!, dvdx_c, v_crop, dx_inv; ndrange = size(dvdx_c), sync = false)
    KA.synchronize(backend)

    launch!(_apply_mask!, dudx, gh.mask; ndrange = size(dudx), sync = false)
    launch!(_apply_mask!, dvdy, gh.mask; ndrange = size(dvdy), sync = false)
    @. shear_c = dudy_c + dvdx_c
    launch!(_apply_mask!, shear_c, gc.mask; ndrange = size(shear_c))
    launch!(_avg_xyT!, shear_h, shear_c; ndrange = size(shear_h))
    launch!(_apply_mask!, shear_h, gh.mask; ndrange = size(shear_h), sync = false)
    KA.synchronize(backend)

    @. gh.shelf_strain_rate = sqrt(dudx^2 + dvdy^2 + dudx * dvdy + 0.25 * shear_h^2)
    return model
end

"""
    update_av_speed!(model::AbstractModel)

Find the depth-averaged speed on the h-grid using components on u- and v- grids
"""
function update_av_speed!(model::AbstractModel)
    @unpack gh,gu,gv=model.fields
    backend = KA.get_backend(gh.h)
    u_crop = similar(gu.u)
    v_crop = similar(gv.v)
    u_h = similar(gh.h)
    v_h = similar(gh.h)

    copyto!(u_crop, gu.u)
    copyto!(v_crop, gv.v)
    launch!(_apply_mask!, u_crop, gu.mask; ndrange = size(u_crop), sync = false)
    launch!(_apply_mask!, v_crop, gv.mask; ndrange = size(v_crop), sync = false)
    KA.synchronize(backend)

    launch!(_avg_x!, u_h, u_crop; ndrange = size(u_h), sync = false)
    launch!(_avg_y!, v_h, v_crop; ndrange = size(v_h), sync = false)
    KA.synchronize(backend)
    launch!(_apply_mask!, u_h, gh.mask; ndrange = size(u_h), sync = false)
    launch!(_apply_mask!, v_h, gh.mask; ndrange = size(v_h), sync = false)
    KA.synchronize(backend)

    @. gh.av_speed = sqrt(u_h^2 + v_h^2)
    return model
end

"""
    update_bed_speed!(model::AbstractModel)

Find the sliding speed at the bed on the h-grid using the average speed.
"""
function update_bed_speed!(model::AbstractModel)
    @unpack gh=model.fields
    gh.bed_speed .= gh.av_speed ./ (1.0 .+ gh.quad_f2 .* gh.β)
    return model
end

"""
    update_β!(model::AbstractModel)

Find the drag coefficient at the bed through the chosen sliding law.
The specific function lives in the corresponding sliding law file.
"""
function update_β!(model::AbstractModel)
    update_β_using_sliding_law!(model.sliding_law,model)
    return model
end


"""
    update_basal_drag!(model::AbstractModel)

Find the shear stress at the bed.
"""
function update_basal_drag!(model::AbstractModel)
    @unpack gh=model.fields
    gh.τbed .= gh.β .* gh.bed_speed
    return model
end



"""
    inner_update_viscosity!(model::AbstractModel)

Inner update to iteratively refine viscosity on the 3d grid at all sigma levels.
"""
@kernel function _inner_update_viscosity_kernel!(
    η,
    glen_b,
    mask,
    shelf_strain_rate,
    τbed,
    ζ,
    glen_reg_strain_rate,
    glen_n_inv_minus_1,
    n_iter_viscosity,
)
    i, j, k = @index(Global, NTuple)
    @inbounds if mask[i, j]
        for iter in 1:n_iter_viscosity
            η[i, j, k] = 0.5 *
                glen_b[i, j, k] *
                (
                    sqrt(
                        shelf_strain_rate[i, j]^2 +
                            0.25 * (τbed[i, j] * ζ[k] / η[i, j, k])^2 +
                            glen_reg_strain_rate^2,
                    )
                )^glen_n_inv_minus_1
        end
    end
end

function inner_update_viscosity!(model::AbstractModel)
    @unpack gh, g3d = model.fields
    @unpack params, solver_params = model
    glen_n_inv_minus_1 = 1.0 / params.glen_n - 1.0
    WAVI.Stencils.launch!(
        _inner_update_viscosity_kernel!,
        g3d.η,
        g3d.glen_b,
        gh.mask,
        gh.shelf_strain_rate,
        gh.τbed,
        g3d.ζ,
        params.glen_reg_strain_rate,
        glen_n_inv_minus_1,
        solver_params.n_iter_viscosity;
        ndrange = (g3d.nxs, g3d.nys, g3d.nσs),
    )
    return model
end

"""
    update_av_viscosity!(model::AbstractModel)

Use quadrature to compute the depth averaged viscosity.
"""
@kernel function _update_av_viscosity_kernel!(ηav, η, mask, quadrature_weights, nσs)
    i, j = @index(Global, NTuple)
    @inbounds if mask[i, j]
        sum_η = zero(eltype(ηav))
        for k in 1:nσs
            sum_η += quadrature_weights[k] * η[i, j, k]
        end
        ηav[i, j] = sum_η
    end
end

function update_av_viscosity!(model::AbstractModel)
    @unpack gh, g3d = model.fields
    gh.ηav .= zero(gh.ηav)
    WAVI.Stencils.launch!(
        _update_av_viscosity_kernel!,
        gh.ηav,
        g3d.η,
        gh.mask,
        g3d.quadrature_weights,
        g3d.nσs;
        ndrange = (g3d.nxs, g3d.nys),
    )
    return model
end

"""
    update_quadrature_falpha!(model::AbstractModel)

Use quadrature to compute falpha functions, used to relate average velocities, basal velocities, and surface velocities to one another
"""
@kernel function _update_quadrature_falpha_kernel!(
    quad_f0,
    quad_f1,
    quad_f2,
    h,
    η,
    ζ,
    mask,
    quadrature_weights,
    nσs,
)
    i, j = @index(Global, NTuple)
    @inbounds if mask[i, j]
        # Creating temp vars for performance (reduced memory traffic)
        f0 = zero(eltype(quad_f0))
        f1 = zero(eltype(quad_f1))
        f2 = zero(eltype(quad_f2))
        h_val = h[i, j]
        for k in 1:nσs
            qw = quadrature_weights[k]
            inv_η = 1.0 / η[i, j, k]
            z_val = ζ[k]
            f0 += qw * h_val * inv_η
            f1 += qw * h_val * z_val * inv_η
            f2 += qw * h_val * (z_val^2) * inv_η
        end
        quad_f0[i, j] = f0
        quad_f1[i, j] = f1
        quad_f2[i, j] = f2
    end
end

function update_quadrature_falpha!(model::AbstractModel)
    @unpack gh, g3d = model.fields
    gh.quad_f0 .= zero(gh.quad_f0)
    gh.quad_f1 .= zero(gh.quad_f1)
    gh.quad_f2 .= zero(gh.quad_f2)
    WAVI.Stencils.launch!(
        _update_quadrature_falpha_kernel!,
        gh.quad_f0,
        gh.quad_f1,
        gh.quad_f2,
        gh.h,
        g3d.η,
        g3d.ζ,
        gh.mask,
        g3d.quadrature_weights,
        g3d.nσs;
        ndrange = (g3d.nxs, g3d.nys),
    )
    return model
end

"""
    update_βeff!(model::AbstractModel)

Compute the effective drag coefficient.
"""
function update_βeff!(model::AbstractModel)
    @unpack gh=model.fields
  #  gh.βeff[gh.mask] .= gh.β[gh.mask] ./ (1.0 .+ gh.quad_f2[gh.mask] .* gh.β[gh.mask])
    gh.βeff .= gh.β ./ (1.0 .+ gh.quad_f2 .* gh.β)
    return model
end



"""
    update_βeff_on_uv_grids!(model::AbstractModel)

Interpolate the effective drag coefficient onto u- and v-grids, accounting for grounded fraction.
"""
function update_βeff_on_uv_grids!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv=model.fields
    @assert eltype(gh.grounded_fraction)==eltype(gh.βeff)
    backend = KA.get_backend(gh.h)

    ones_crop = similar(gh.h)
    β_crop = similar(gh.h)
    gf_crop = similar(gh.h)
    tmpu = similar(gu.u)
    tmpv = similar(gv.v)
    denu = similar(gu.u)
    denv = similar(gv.v)
    ipolgfu=zeros(T,gu.nxu,gu.nyu)
    ipolgfv=zeros(T,gv.nxv,gv.nyv)

    fill!(ones_crop, one(T))
    copyto!(β_crop, gh.βeff)
    copyto!(gf_crop, gh.grounded_fraction)
    launch!(_apply_mask!, ones_crop, gh.mask; ndrange = size(ones_crop), sync = false)
    launch!(_apply_mask!, β_crop, gh.mask; ndrange = size(β_crop), sync = false)
    launch!(_apply_mask!, gf_crop, gh.mask; ndrange = size(gf_crop), sync = false)
    KA.synchronize(backend)

    launch!(_avg_xT!, denu, ones_crop; ndrange = size(denu), sync = false)
    launch!(_avg_yT!, denv, ones_crop; ndrange = size(denv), sync = false)
    KA.synchronize(backend)

    launch!(_avg_xT!, tmpu, β_crop; ndrange = size(tmpu))
    @views gu.βeff[gu.mask] .= tmpu[gu.mask] ./ denu[gu.mask]
    launch!(_avg_xT!, tmpu, gf_crop; ndrange = size(tmpu))
    @views ipolgfu[gu.mask] .= tmpu[gu.mask] ./ denu[gu.mask]
    gu.βeff[ipolgfu .> zero(T)] .= gu.βeff[ipolgfu .> zero(T)].*gu.grounded_fraction[ipolgfu .> zero(T)]./
                                                        ipolgfu[ipolgfu .> zero(T)]

    launch!(_avg_yT!, tmpv, β_crop; ndrange = size(tmpv))
    @views gv.βeff[gv.mask] .= tmpv[gv.mask] ./ denv[gv.mask]
    launch!(_avg_yT!, tmpv, gf_crop; ndrange = size(tmpv))
    @views ipolgfv[gv.mask] .= tmpv[gv.mask] ./ denv[gv.mask]
    gv.βeff[ipolgfv .> zero(T)] .= gv.βeff[ipolgfv .> zero(T)].*gv.grounded_fraction[ipolgfv .> zero(T)]./
                                                 ipolgfv[ipolgfv .> zero(T)];

    return model
end


"""
    update_rheological_operators!(model::AbstractModel)

Precompute various diagonal matrices used in defining the momentum operator.
"""
function update_rheological_operators!(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc = model.fields
    @unpack params, solver_params = model

    hη = similar(gh.h)
    hη_c = similar(gh.h, T, gc.nxc, gc.nyc)
    @. hη = gh.h * gh.ηav
    launch!(_avg_xy!, hη_c, hη; ndrange = size(hη_c))

    gh_diag = reshape(gh.dneghηav[].diag, gh.nxh, gh.nyh)
    gc_diag = reshape(gc.dneghηav[].diag, gc.nxc, gc.nyc)
    gu_diag = reshape(gu.dnegβeff[].diag, gu.nxu, gu.nyu)
    gv_diag = reshape(gv.dnegβeff[].diag, gv.nxv, gv.nyv)
    gh_imp = reshape(gh.dimplicit[].diag, gh.nxh, gh.nyh)

    @. gh_diag = ifelse(gh.mask, -gh.h * gh.ηav, zero(T))
    @. gc_diag = ifelse(gc.mask, -hη_c, zero(T))
    @. gu_diag = ifelse(gu.mask, -gu.βeff, zero(T))
    @. gv_diag = ifelse(gv.mask, -gv.βeff, zero(T))
    @. gh_imp = ifelse(
        gh.mask,
        -params.density_ice * params.g * solver_params.super_implicitness * params.dt * gh.dsdh,
        zero(T),
    )
    return model
end




"""
    get_op(model::AbstractModel{T,N}) where {T,N}

 Get operator, defined as a LinearMap type.

"""
function get_op(model::AbstractModel{T,N}) where {T,N}
    @unpack gu,gv=model.fields
    ni = gu.ni + gv.ni
    op_fun! = get_op_fun(model)
    op=LinearMap{T}(op_fun!,ni;issymmetric=true,ismutating=true,ishermitian=true,isposdef=true)
end


"""
    get_rhs_dirichlet(model::AbstractModel{T,N}) where {T,N}

    Extra term of right hand side to implement non-homogenous Dirichlet conditions   

"""
function get_rhs_dirichlet!(rhs_dirichlet,model::AbstractModel{T,N}) where {T,N}
    @unpack gu,gv=model.fields

    uvfixed=[
    gu.u[:].*gu.u_isfixed[:]
    ;
    gv.v[:].*gv.v_isfixed[:]
    ]

    op_fun! = get_op_fun(model)
    op_fun!(rhs_dirichlet,uvfixed,vecSampled=false)
    
    @. rhs_dirichlet = - rhs_dirichlet
    
    return rhs_dirichlet
end

"""
    set_residual!(model::AbstractModel,residual)

Set residuals to particular values. Input vector residual represents stacked u and v components at valid grid points.
"""
function set_residual!(model::AbstractModel,residual)
    @unpack gu,gv=model.fields
    @views gu.residual[gu.mask_inner] .= residual[1:gu.ni]
    @views gv.residual[gv.mask_inner] .= residual[(gu.ni+1):(gu.ni+gv.ni)]
    return model
end
