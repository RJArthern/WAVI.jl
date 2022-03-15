using WAVI 
function MISMIP_PLUS()
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
    bed = WAVI.mismip_plus_bed #function definition

    #solver parameters
    maxiter_picard = 1
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    default_thickness = 100.0 #set the initial condition this way
    accumulation_rate = 0.3
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate)

    #make the model
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     solver_params = solver_params)

    #timestepping parameters
    dt = 1.0
    end_time = 1.0
    timestepping_params = TimesteppingParams(niter0 = 0, 
                                            dt = dt, 
                                            end_time = end_time);

   
    
    @time simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params);
            
    #perform the simulation
    @time run_simulation!(simulation);
    return simulation;
end

t = @elapsed MISMIP_PLUS(); println("total time: ", t)