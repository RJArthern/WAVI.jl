"""
Iceberg(n_timesteps)

Driver routine for iceberg simulation using WAVI.
Defines parameters, initialises state, and runs n_timesteps.
"""
function iceberg_test(;end_time=1000.)

#Number of timesteps between plots
n_steps_plot=Inf

# Iceberg
nx = 20
ny = 20
nσ = 4
x0 = -50000.0
y0 = -50000.0
dx = 5000.0
dy = 5000.0
#Homogenous Dirichlet boundary conditions
u_iszero=falses(nx+1,ny)
u_iszero[div(nx,2)+1,:].=true
v_iszero=falses(nx,ny+1)
v_iszero[:,div(ny,2)+1].=true
grid = Grid(nx = nx, 
            ny = ny, 
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            u_iszero = u_iszero,
            v_iszero = v_iszero)

h_mask=sqrt.(grid.xxh.^2+grid.yyh.^2) .< 45000
h_mask = convert(Array{Bool,2}, h_mask)
grid = @set grid.h_mask = h_mask


#Model domain mask on h-grid
starting_thickness=zeros(nx,ny)
starting_thickness[h_mask] .= 200.0
initial_conditions = InitialConditions(initial_thickness = starting_thickness)

accumulation_rate=0.3
params = Params(accumulation_rate = accumulation_rate)

bed_elevation=-500.0.*ones(nx,ny)

maxiter_picard=1 #No need for Picard iteration if runninng to steady state
solver_params = SolverParams(maxiter_picard = maxiter_picard)

#make the model
model = Model(grid = grid,
              bed_elevation = bed_elevation, 
              params = params, 
              solver_params = solver_params, 
              initial_conditions = initial_conditions)

timestepping_params = TimesteppingParams(dt = 0.1, end_time = end_time)
simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params)
run_simulation!(simulation)
return simulation
end
