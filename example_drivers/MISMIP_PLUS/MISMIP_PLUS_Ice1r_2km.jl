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

    #Initial conditions from Ice0 
    h_init = Array{Float64}(undef, nx, ny);
    read!("examples\\MISMIP_ice0_final.bin", h_init)
    h_init = ntoh.(h_init);
    initial_conditions = InitialConditions(initial_thickness = h_init)

    #solver parameters
    maxiter_picard = 30
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    accumulation_rate = 0.3
    default_temperature=265.700709
    params = Params(accumulation_rate = accumulation_rate,
                    default_temperature = default_temperature)

    #make the model
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     solver_params = solver_params, 
                     initial_conditions = initial_conditions)

    #embed the ice model with melt rate model
    function m1(h, b)
        draft = -(918.0 / 1028.0) * h
        cavity_thickness = draft .- b
        cavity_thickness = max.(cavity_thickness, 0)
        m =  0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
        return m
    end
    arguments = (h = model.fields.gh.h, b = model.fields.gh.b)
    melt_rate_model = AnalyticMeltRate(melt_partial_cell = true, 
                                    melt_rate_function = m1, 
                                    function_arguments = arguments)
    add_melt_rate_model!(model,melt_rate_model)

    #timestepping parameters
    niter0 = 0
    dt = 0.001
    end_time = 1.
    chkpt_freq = 0.1
    pchkpt_freq = 0.1
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
    output_freq = 0.1
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

@time simulation = MISMIP_PLUS_Ice1r_2km();