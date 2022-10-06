
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
tidal_lengthscale=5*1e3
tidal_melting = false
tidal_drag = true
default_thickness = 680.0 #set the initial condition this way
accumulation_rate = 0.3
params = Params(default_thickness = default_thickness, 
               accumulation_rate = accumulation_rate,
               tidal_drag = tidal_drag,
               tidal_lengthscale = tidal_lengthscale,
               tidal_melting = tidal_melting)

#make the model
model = Model(grid = grid,
                    bed_elevation = bed, 
                    params = params, 
                    solver_params = solver_params)

                
dt = 0.5;
end_time = 0.5;
timestepping_params = TimesteppingParams(dt = dt, 
                                         end_time = end_time)

simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params)          
                
run_simulation!(simulation)

#plot the basal friction field
heatmap(simulation.model.grid.xxh[:,1],simulation.model.grid.yyh[1,:],model.fields.gh.weertman_c')
contour!(simulation.model.grid.xxh[:,1],simulation.model.grid.yyh[1,:], model.fields.gh.grounded_fraction', levels=[0.5,0.5], linecolor = :blue)