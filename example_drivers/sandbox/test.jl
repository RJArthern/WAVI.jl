
#Grid and boundary conditions
nx = 160
ny = 20
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 4000.0
dy = 4000.0
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
#
##Physical parameters
tidal_lengthscale=4*1e3
tidal_melting = true
tidal_drag = true
partial_cell_drag = false #only relevant if tidal drag = false
tidal_daily_timescale = 14.0; 
tidal_hourly_timescale = 6.25; 
#to avoid any oscillations
tidal_daily_timescale = 1.0e35; 
tidal_hourly_timescale = 1.0e35; 

default_thickness = 680.0 #set the initial condition this way
accumulation_rate = 0.3
params = Params(default_thickness = default_thickness, 
               accumulation_rate = accumulation_rate,
               tidal_drag = tidal_drag,
               tidal_lengthscale = tidal_lengthscale,
               tidal_melting = tidal_melting,
               partial_cell_drag = partial_cell_drag, 
               tidal_daily_timescale = tidal_daily_timescale, 
               tidal_hourly_timescale = tidal_hourly_timescale)

#make the model
melt_rate = UniformMeltRate(m = 30, partial_cell_melting = false)
model = Model(grid = grid,
                    bed_elevation = bed, 
                    params = params, 
                    solver_params = solver_params,
                    melt_rate = melt_rate)

                
dt = 1/(365*24);
end_time = 20*dt
timestepping_params = TimesteppingParams(dt = dt, 
                                         end_time = end_time)

simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params)          
                
run_simulation!(simulation)

#plot the basal friction or melt field
heatmap(simulation.model.grid.xxh[:,1],simulation.model.grid.yyh[1,:],model.fields.gh.weertman_c')
#heatmap(simulation.model.grid.xxh[:,1],simulation.model.grid.yyh[1,:],model.fields.gh.basal_melt')
contour!(simulation.model.grid.xxh[:,1],simulation.model.grid.yyh[1,:], model.fields.gh.grounded_fraction', levels=[0.5,0.5], linecolor = :blue)