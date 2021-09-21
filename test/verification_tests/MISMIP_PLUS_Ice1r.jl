using WAVI 
function MISMIP_PLUS_Ice1r(initial_thickness)
    #Grid and boundary conditions
    nx = 80
    ny = 10
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 8000.0
    dy = 8000.0
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
    bed = WAVI.mismip_plus_bed #function definition

    #solver parameters
    maxiter_picard = 30
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    accumulation_rate = 0.3
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate)

    #initial conditions
    initial_conditions = InitialConditions(initial_thickness = initial_thickness)

    #make the model
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     solver_params = solver_params,
                     initial_conditions = initial_conditions)

    #embed the model with melt rate model
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
    dt = 0.1
    end_time = 100.
    timestepping_params = TimesteppingParams(dt = dt, end_time = end_time)
    
    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

#@time simulation = MISMIP_PLUS_test();