using WAVI, Plots
function plot_MISMIP1_melt_pattern()
    # Grid and boundary conditions
    nx = 320
    ny = 40
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 2000.0
    dy = 2000.0
    h_mask = trues(nx, ny)
    u_iszero = falses(nx + 1, ny); u_iszero[1,:] .= true
    v_iszero = falses(nx, ny + 1); v_iszero[:,1] .= true; v_iszero[:,end] .= true
    grid = Grid(nx=nx, 
                ny=ny,   
                nσ=nσ, 
                x0=x0, 
                y0=y0, 
                dx=dx, 
                dy=dy,
                h_mask=h_mask, 
                u_iszero=u_iszero, 
                v_iszero=v_iszero)

    # Bed 
    bed = WAVI.mismip_plus_bed # function definition

    # Inputing thickness profile, reading from binary file 
    fname = "examples\\Favier2019_melt_params\\data\\MISMIP_ice0_2km_SteadyThickness.bin";
    # fname = joinpath(dirname(@__FILE__), "data",  "MISMIP_ice0_2km_SteadyThickness.bin")
    h = Array{Float64,2}(undef, nx, ny)
    read!(fname, h)
    h = ntoh.(h)

    # melt rate Model    
    melt_rate = MISMIPMeltRateOne();

    # make the model
    initial_conditions = InitialConditions(initial_thickness=h) # set thickness 
    model = Model(grid=grid,
                bed_elevation=bed, 
                initial_conditions=initial_conditions,
                melt_rate=melt_rate,
                solver_params=SolverParams(maxiter_picard=1) # we'll update state to get basal melt, but don't need other info accurate
                ) 

    # update the configuration
    update_state!(model)
    
    # make the plot
    x = model.grid.xxh[:,1]; y = model.grid.yyh[1,:];
    m = deepcopy(model.fields.gh.basal_melt);
    m[model.fields.gh.grounded_fraction .== 1.] .= NaN;
    msat = deepcopy(m)
    msat[msat .> 50] .= 50
    plt = contour(x ./ 1e3,y / 1e3, msat', 
                  fill=true, 
                  linewidth=0, 
                  colorbar=true,
                  colorbar_title="melt rate (m/yr)",
                  framestyle=:box);
    xlims!((420, 600))
    xlabel!("x (km)")
    ylabel!("y (km)")
    @show plt

    #savefig(plt,  joinpath(dirname(@__FILE__), "MISMIP1_melt_pattern.png"))

    mean_melt_rate = sum(m[model.fields.gh.grounded_fraction .< 1]) / sum(model.fields.gh.grounded_fraction .< 1)
    println("Mean melt rate in the shelf is ", mean_melt_rate, " m/yr")
  
    return model
end

model = plot_MISMIP1_melt_pattern()
