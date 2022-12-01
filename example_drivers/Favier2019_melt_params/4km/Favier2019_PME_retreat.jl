using WAVI 
function MISMIP_PLUS_Ice1r_2km()
    #Grid and boundary conditions
    nx = 320
    ny = 40
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 2000.0
    dy = 2000.0
    h_mask=trues(nx,ny)
    u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
    v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true
    grid = Grid(nx = nx, 
                ny = ny,   
                nσ = nσ, 
                x0 = x0, 
                y0 = y0, 
                dx = dx, 
                dy = dy,
                h_mask = h_mask, 
                u_iszero = u_iszero, 
                v_iszero = v_iszero)

    #Bed 
    bed = WAVI.mismip_plus_bed #function definition of bed

    # Inputing thickness profile, reading from binary file 
    fname = "examples\\Favier2019_melt_params\\data\\MISMIP_ice0_2km_SteadyThickness.bin";
    # fname = joinpath(dirname(@__FILE__), "data",  "MISMIP_ice0_2km_SteadyThickness.bin")
    h = Array{Float64,2}(undef, nx, ny)
    read!(fname, h)
    h = ntoh.(h)

    #solver parameters
    maxiter_picard = 30
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    # melt rate Model    
    pme = PlumeEmulator(α=1.5)

    #Physical parameters
    accumulation_rate = 0.3
    params = Params(accumulation_rate = accumulation_rate)

    #make the model
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     melt_rate = pme,
                     solver_params = solver_params, 
                     initial_conditions = initial_conditions)



    #timestepping parameters
    niter0 = 0
    dt = 0.02
    end_time = 100.
    chkpt_freq = 1.
    pchkpt_freq = 1.
    timestepping_params = TimesteppingParams(niter0 = niter0, 
                                            dt = dt, 
                                            end_time = end_time, 
                                            chkpt_freq = chkpt_freq, 
                                            pchkpt_freq = pchkpt_freq)

    #output parameters
    folder = "outputs"
    isdir(folder) && rm(folder, force = true, recursive = true)
    mkdir(folder) #make a clean folder for outputs
    outputs = (h   = model.fields.gh.h,
                u  = model.fields.gh.u,
                v  = model.fields.gh.v,
                melt = model.fields.gh.basal_melt,
                grfrac = model.fields.gh.grounded_fraction)
    output_freq = 1.
    output_params = OutputParams(outputs = outputs, 
                            output_freq = output_freq,
                            output_format = "mat",
                            output_path = folder,
                            zip_format = "nc")
    
    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

@time simulation = Favier2019_PME_retreat();