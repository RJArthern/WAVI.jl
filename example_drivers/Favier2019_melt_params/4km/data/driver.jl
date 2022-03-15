using WAVI

function driver()
# Grid and boundary conditions
nx = 160
ny = 20
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 4000.0
dy = 4000.0
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
fname = ("WAVI_ice0_4km_thickness.bin")
h = Array{Float64,2}(undef, nx, ny)
read!(fname, h)
h = ntoh.(h)

# make the model
initial_conditions = InitialConditions(initial_thickness=h) # set thickness

#specify melt model
 melt_model = QuadraticMeltRate(γT = 0.4*1e-3)

model = Model(grid=grid,
            bed_elevation=bed,
            melt_rate = melt_model,
            initial_conditions=initial_conditions,
            solver_params=SolverParams(maxiter_picard=10) ,
            params=Params(accumulation_rate=0.3)
            ) ;


#timestepping parameters
niter0 = 0
dt = 0.02
end_time = 100.
chkpt_freq = 5.
pchkpt_freq = 5.
timestepping_params = TimesteppingParams(niter0 = niter0,
                                        dt = dt,
                                        end_time = end_time,
                                        chkpt_freq = chkpt_freq,
                                        pchkpt_freq = pchkpt_freq)



#output parameters
outputs = (h   = model.fields.gh.h,
            u  = model.fields.gh.u,
            v  = model.fields.gh.v,
            melt = model.fields.gh.basal_melt,
            grfrac = model.fields.gh.grounded_fraction)
output_freq = 0.01
output_params = OutputParams(outputs = outputs,
                        output_freq = output_freq,
                        output_format = "jld2",
                        zip_format = "nc",
                        output_start = true)



                        simulation = Simulation(model = model,
                        timestepping_params = timestepping_params,
                        output_params = output_params);
    
    #perform the simulation
    run_simulation!(simulation);
    
    return simulation
    end
    
    