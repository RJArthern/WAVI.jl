# # MISMIP+ example
# 
# This example demonstrates the Marine Ice Sheet Model Intercomparison (MISMIP) + ice0 experiment (doi: 10.5194/tc-14-2283-2020). 
# This intercomparison exercise considers a plan-view ice sheet, in which the grounding line can stabilize on a section of bed which has a locally positive slope in the flow direction. 
# Such configurations are theoretically impossible in one horizontal dimension (see doi: 10.1029/2006JF000664), demonstrating the importance of variations in the second dimension for buttressing ice sheets. 
# This example demonstrates how to apply boundary conditions in `WAVI.jl`, and control the number of iterations in the velocity solve. 

# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add WAVI, Plots"
# ```
using WAVI, Plots

# ## Basal Topography
#
# The MISMIP+ domain is 640km in the x-direction and 80km in the y-direction, centred around $y = 0$.
# The basal topography is given by $z_b = \max [B_x(x) + B_y(y), -720] where $B_x(x)$ is a sixth order, even polynomial and $B_y(y)$ introduces two bumps in the domain:
function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75
    wc = 24000.0; fc = 4000.0; dc = 500.0
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) )
    b = max(bx(x) + by(y), -720.0)
    return b
end

# Let's take a look at this bed. First we define the grid sizes and build some arrays, so we can plot. We'll use a high resolution to get a nice plot:
dx = 1.e3
dy = 1.e3
nx = round(Int, 640*1e3/dx)
ny = round(Int, 80*1e3/dx)
xx=[i*dx for i=1:nx, j=1:ny]
yy=[j*dy for i=1:nx, j=1:ny] .- 42000
x = xx[:,1]
y = yy[1,:]

# Now we can plot
plt =  Plots.heatmap(x/1e3, y/1e3, mismip_plus_bed.(xx,yy)', 
                    xlabel = "x (km)", 
                    ylabel = "y (km)",
                    colorbar_title = "bed depth (m)")
plot!(size = (800,400))

# ## Boundary Conditions 
#
# The MISMIP+ experiment specifies no slip (zero velocity in both directions) boundary conditions at $x = 0$, and free-slip boundary conditions (zero velocity in the direction normal to the walls) at the lateral boundaries at $y = 0$km and $y = 84$km.
# Velocity boundary conditions are controlled by specifying zeros in appropriate entries in arrays, which are then passed to the grid:
u_iszero = falses(nx+1,ny); #build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere 
u_iszero[1,:].=true;        #set the x-direction velocity to zero at x = 0.
v_iszero=falses(nx,ny+1);   #build x-direction velocity boundary condition matrix with no zero boundary conditions anywhere 
v_iszero[:,1].=true;        #set the y-direction velocity to zero at y = 0 (free slip)
v_iszero[:,end].=true       #set the y-direction velocity to zero at y = 84km (free slip)
v_iszero[1,:].=true         #set the y-direction velocity to zero at x = 0km (no slip in combination with u_iszero)

# Now we build the grid as usual, passing the arrays we just constructed via optional arguments. 
# We'll also redefine the grip size to be lower resolution to make the later simulations quicker
dx = 4.e3
dy = 4.e3
nx = round(Int, 640*1e3/dx)
ny = round(Int, 80*1e3/dx)
grid = Grid(nx = nx, 
            ny = ny,   
            dx = dx, 
            dy = dy,
            u_iszero = u_iszero, 
            v_iszero = v_iszero)

# ## Solver Parameters
# We're interested in the steady state reached at late times, rather than the solution along the way. We don't need to get the velocity right along the way, just have it correct eventually.
# We therefore set the number of iterations in the velocity solve to be small: at each timestep, the solver just does a small number of iterations, and the velocity is only approximate. But, since we do a lot of iterations getting to steady state, the velocity gets to the right thing eventually.
# This number, and other parameters relating to how the equations are solved, are set via a `SolverParams` object:
solver_params = SolverParams(maxiter_picard = 1)
# Explicitly, we set the number of iterations (formally, Picard iterations) to be as small as possible, i.e. one iteration.

# ## Make the model
# Now we have our grid, bed, and solver parameters, we just need to set the appropriate initial conditions and physical parameters for MISMIP+, and then we can build our model. 
# In MISMIP+, the initial condition specifies 100m thick ice everywhere
initial_conditions = InitialConditions(initial_thickness = 100 .* ones(nx,ny))

# And the accumulation rate is set to 0.3m/a (all other defaults in WAVI are chosen according to the values in MISMIP+)
params = Params( accumulation_rate = accumulation_rate)

# Now let's make our model! Note that we use the functional form of the bed (the array we plotted earlier has higher resolution than our model has)
model = Model(grid = grid,
            bed_elevation = mismip_plus_bed, 
            initial_conditions = initial_conditions,
            params = params, 
            solver_params = solver_params)


