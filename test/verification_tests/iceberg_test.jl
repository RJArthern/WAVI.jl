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

#Inhomogenous Dirichlet boundary conditions will be applied on x- and y- center lines
u_isfixed=falses(nx+1,ny)
u_isfixed[div(nx,2)+1,:].=true
v_isfixed=falses(nx,ny+1)
v_isfixed[:,div(ny,2)+1].=true

grid = Grid(nx = nx, 
            ny = ny, 
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            u_isfixed = u_isfixed,
            v_isfixed = v_isfixed)

h_mask=sqrt.(grid.xxh.^2+grid.yyh.^2) .< 45000
h_mask = convert(Array{Bool,2}, h_mask)
grid = @set grid.h_mask = h_mask


#Model domain mask on h-grid
starting_thickness=zeros(nx,ny)
starting_thickness[h_mask] .= 200.0

u_drift=5.0
v_drift=3.0
rotation_period_yrs = 10000.0
rotation_rate_per_yr = 2*pi/rotation_period_yrs 
solid_body_u = u_drift .* ones(nx+1,ny) .- rotation_rate_per_yr .* grid.yyu
solid_body_v = v_drift .* ones(nx,ny+1) .+ rotation_rate_per_yr .* grid.xxv


u_bc = zeros(nx+1,ny)
v_bc = zeros(nx,ny+1)

u_bc[u_isfixed] .= solid_body_u[u_isfixed]
v_bc[v_isfixed] .= solid_body_v[v_isfixed]

initial_conditions = InitialConditions(initial_thickness = starting_thickness, 
                                       initial_u_veloc = u_bc,
                                       initial_v_veloc = v_bc)

accumulation_rate=0.3
default_temperature=265.700709
params = Params(accumulation_rate = accumulation_rate, default_temperature = default_temperature)

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

h0=((36.0*0.3/(2.0e-17))*(1.0/(9.81*918.0*(1-918.0/1028.0)))^3)^(1.0/4.0)
u0 = simulation.model.grid.xxu*0.3/(2.0*h0) .+ solid_body_u
v0 = simulation.model.grid.yyv*0.3/(2.0*h0) .+ solid_body_v
relerr_h=norm(simulation.model.fields.gh.h[simulation.model.fields.gh.mask].-h0)/
            norm(h0*ones(length(simulation.model.fields.gh.h[simulation.model.fields.gh.mask])))
relerr_u=norm(simulation.model.fields.gu.u[simulation.model.fields.gu.mask]-u0[simulation.model.fields.gu.mask])/
                 norm(u0[simulation.model.fields.gu.mask])
relerr_v=norm(simulation.model.fields.gv.v[simulation.model.fields.gv.mask]-v0[simulation.model.fields.gv.mask])/
                 norm(v0[simulation.model.fields.gv.mask])

return simulation, relerr_h, relerr_u, relerr_v

end
