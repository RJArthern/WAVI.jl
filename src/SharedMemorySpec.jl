

@with_kw struct SharedMemorySpec{T,N} <: AbstractParallelSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    schwarzModelArray::Array{Model,2} = Array{Model,2}(undef,ngridsx,ngridsy)
end

get_parallel_spec(model::AbstractModel) = model.parallel_spec

update_preconditioner!(model::AbstractModel) = update_preconditioner!(model::AbstractModel,get_parallel_spec(model::AbstractModel))

function update_preconditioner!(model::AbstractModel,::BasicParallelSpec)
    return model
end

function update_preconditioner!(model::AbstractModel,::SharedMemorySpec)
    @unpack ngridsx, ngridsy, overlap = model.parallel_spec

    @sync for igrid = 1:ngridsx
        for jgrid = 1:ngridsy
            Threads.@spawn begin
                model.parallel_spec.schwarzModelArray[igrid,jgrid] = schwarzModel(model;
                                                                                  igrid=igrid,
                                                                                  jgrid=jgrid,
                                                                                  ngridsx=ngridsx,
                                                                                  ngridsy=ngridsy, 
                                                                                  overlap=overlap)
            end
        end
    end
    return model
end


function precondition!(model::AbstractModel,::SharedMemorySpec)

@unpack ngridsx, ngridsy, overlap, niterations, schwarzModelArray, damping = model.parallel_spec
@unpack solver_params = model

x=get_start_guess(model)
    
op=get_op(model)

b=get_rhs(model)

resid=get_resid(x,op,b)

set_residual!(model,resid)

rel_resid = norm(resid)/norm(b)

converged = rel_resid < solver_params.tol_picard

if ! converged

    for iteration = 1:niterations
        println("")
        @sync for igrid = 1:ngridsx
            for jgrid = 1:ngridsy
                Threads.@spawn begin
            
                    model_g = schwarzModelArray[igrid,jgrid]

                    schwarzRestrictVelocities!(model_g::AbstractModel,
                                            model::AbstractModel;
                                            igrid=igrid,
                                            jgrid=jgrid,
                                            ngridsx=ngridsx,
                                            ngridsy=ngridsy,
                                            overlap=overlap)
                end
            end
        end
        
        @sync for igrid = 1:ngridsx
            for jgrid = 1:ngridsy
                Threads.@spawn begin
                    model_g = schwarzModelArray[igrid,jgrid]
                    update_state!(model_g)
                end
            end
        end

        model.fields.gu.u[:,:] .= damping .* model.fields.gu.u
        model.fields.gv.v[:,:] .= damping .* model.fields.gv.v
        
        threadLock=ReentrantLock()
        @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                        Threads.@spawn begin

                                model_g = schwarzModelArray[igrid,jgrid]

                                lock(threadLock)
                                try
                                    schwarzProlongVelocities!(model::AbstractModel,
                                                            model_g::AbstractModel;
                                                            igrid=igrid,
                                                            jgrid=jgrid,
                                                            ngridsx=ngridsx,
                                                            ngridsy=ngridsy,
                                                            overlap=overlap,
                                                            damping=damping)
                                finally
                                    unlock(threadLock)
                                end
                        end                    
                end
        end
        println("")
    end

end
return converged, rel_resid
end


function schwarzModel(model::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1)
    @unpack nx,ny,dx,dy,nσ,x0,y0,h_mask,h_isfixed,u_iszero,v_iszero,u_isfixed,v_isfixed,quadrature_weights,σ = model.grid
    @unpack gh,gu,gv,g3d = model.fields
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    nx_g = i_stop_g - i_start_g + 1
    ny_g = j_stop_g - j_start_g + 1
    dx_g = dx
    dy_g = dy
    nσ_g = nσ
    x0_g = x0 + (i_start_g-1)*dx
    y0_g = y0 + (j_start_g-1)*dy
    h_mask_g = h_mask[i_start_g:i_stop_g,j_start_g:j_stop_g]
    h_isfixed_g = h_isfixed[i_start_g:i_stop_g,j_start_g:j_stop_g]

    u_iszero_g = u_iszero[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    v_iszero_g = v_iszero[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    u_isfixed_g = u_isfixed[i_start_g:i_stop_g+1,j_start_g:j_stop_g] 
    v_isfixed_g = v_isfixed[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    #Set halo points as fixed velocity points
    (igrid==1) || (u_isfixed_g[1,:] .= true)
    (igrid==ngridsx) || (u_isfixed_g[nx_g+1,:] .= true)
    (igrid==1) || (v_isfixed_g[1,:] .= true)
    (igrid==ngridsx) || (v_isfixed_g[nx_g,:] .= true)
    (jgrid==1) || (v_isfixed_g[:,1] .= true)
    (jgrid==ngridsy) || (v_isfixed_g[:,ny_g+1] .= true)
    (jgrid==1) || (u_isfixed_g[:,1] .= true)
    (jgrid==ngridsy) || (u_isfixed_g[:,ny_g] .= true)

    quadrature_weights_g = quadrature_weights 
    σ_g = σ 

    grid_g=Grid(
    nx=nx_g,
    ny=ny_g,
    dx=dx_g,
    dy=dy_g,
    nσ=nσ_g,
    x0=x0_g,
    y0=y0_g,
    h_mask = h_mask_g,
    h_isfixed = h_isfixed_g,
    u_iszero = u_iszero_g,
    v_iszero = v_iszero_g,
    u_isfixed = u_isfixed_g,
    v_isfixed = v_isfixed_g,
    quadrature_weights = quadrature_weights_g,
    σ = σ_g)

    bed_elevation_g = model.fields.gh.b[i_start_g:i_stop_g,j_start_g:j_stop_g]

    params_g = model.params
    params_g = @set params_g.weertman_c = params_g.weertman_c[i_start_g:i_stop_g,j_start_g:j_stop_g]
    params_g = @set params_g.accumulation_rate = params_g.accumulation_rate[i_start_g:i_stop_g,j_start_g:j_stop_g]
    params_g = @set params_g.glen_a_ref = params_g.glen_a_ref[i_start_g:i_stop_g,j_start_g:j_stop_g]

    solver_params_g=model.solver_params

    initial_thickness_g = gh.h[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_grounded_fraction_g = gh.grounded_fraction[i_start_g:i_stop_g,j_start_g:j_stop_g]
    initial_u_veloc_g = gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    initial_v_veloc_g = gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1]
    initial_viscosity_g = g3d.η[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_temperature_g = g3d.θ[i_start_g:i_stop_g,j_start_g:j_stop_g,:]
    initial_damage_g = g3d.Φ[i_start_g:i_stop_g,j_start_g:j_stop_g,:]

    initial_conditions_g=InitialConditions(
        initial_thickness = initial_thickness_g,
        initial_grounded_fraction = initial_grounded_fraction_g,
        initial_u_veloc = initial_u_veloc_g,
        initial_v_veloc = initial_v_veloc_g,
        initial_viscosity = initial_viscosity_g,
        initial_temperature = initial_temperature_g,
        initial_damage = initial_damage_g)

    melt_rate_g=model.melt_rate

    parallel_spec_g = BasicParallelSpec()

    model_g=Model(
        grid = grid_g, 
        bed_elevation = bed_elevation_g,
        params = params_g,
        solver_params = solver_params_g,
        initial_conditions = initial_conditions_g,
        melt_rate = melt_rate_g,
        parallel_spec = parallel_spec_g)

    return model_g
end


function schwarzRestrictVelocities!(model_g::AbstractModel,model::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1)
    @unpack nx,ny = model.grid
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    model_g.fields.gu.u .= model.fields.gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g]
    model_g.fields.gv.v .= model.fields.gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1]

    return model_g
end

function schwarzProlongVelocities!(model::AbstractModel,model_g::AbstractModel;igrid=1,jgrid=1,ngridsx=1,ngridsy=1,overlap=1,damping=0.0)
    @unpack nx,ny = model.grid
    
    @assert rem(nx,ngridsx)==0 "Model domain is not an integer number of subdomains in x-direction"
    @assert rem(ny,ngridsy)==0 "Model domain is not an integer number of subdomains in y-direction"
    (1 <= igrid <= ngridsx) || throw(BoundsError(model,igrid))
    (1 <= jgrid <= ngridsy) || throw(BoundsError(model,jgrid))

    nx_domain=div(nx,ngridsx)
    ny_domain=div(ny,ngridsy)
    i_start_domain=(igrid-1)*nx_domain+1
    i_stop_domain=igrid*nx_domain
    j_start_domain=(jgrid-1)*ny_domain+1
    j_stop_domain=jgrid*ny_domain

    i_start_g = max(i_start_domain - overlap, 1)
    i_stop_g = min(i_stop_domain + overlap, nx)
    j_start_g = max(j_start_domain - overlap, 1)
    j_stop_g = min(j_stop_domain + overlap, ny)

    nx_g = i_stop_g - i_start_g + 1
    ny_g = j_stop_g - j_start_g + 1

    upou = schwarzPartitionOfUnity(nx_g+1,ny_g,igrid==1,igrid==ngridsx,jgrid==1,jgrid==ngridsy,2*overlap,2*overlap-1)
    vpou = schwarzPartitionOfUnity(nx_g,ny_g+1,igrid==1,igrid==ngridsx,jgrid==1,jgrid==ngridsy,2*overlap-1,2*overlap)

    model.fields.gu.u[i_start_g:i_stop_g+1,j_start_g:j_stop_g] .+= (one(damping)-damping) .* upou .* model_g.fields.gu.u
    model.fields.gv.v[i_start_g:i_stop_g,j_start_g:j_stop_g+1] .+= (one(damping)-damping) .* vpou .* model_g.fields.gv.v

    return model
end

function schwarzPartitionOfUnity(m,n,leavei1,leaveim,leavej1,leavejn,overlapi,overlapj)
    @assert m > (~leavei1 && overlapi) + (~leaveim && overlapi)
    @assert n > (~leavej1 && overlapj) + (~leavejn && overlapj)
    pou = [min(1.0, 
        leavei1 ? Inf : (i-1)./(overlapi),
        leaveim ? Inf : (m-i)./(overlapi)) for i=1:m, j=1:n] .*
        [min(1.0, 
        leavej1 ? Inf : (j-1)./(overlapj),
        leavejn ? Inf : (n-j)./(overlapj)) for i=1:m, j=1:n]
    return pou
end