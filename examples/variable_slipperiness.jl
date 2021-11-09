# # Ice Flow over a bed with variable slipperiness and viscosity.

# In this example, we consider flow in two horizontal dimensions, with variable ice viscosity and basal slipperiness. 
# This example is similar in spirit to the MISMIP 3D Intercomparison exercise (http://homepages.ulb.ac.be/~fpattyn/mismip3d/)
# This example demonstrates how to apply spatially heterogenous basal slipperiness and ice viscosity fields.

# ## Install dependencies
#
# First let's make sure we have all required packages installed. 

#using Pkg
#Pkg.add("WAVI"); Pkg.add("Plots"); 
using WAVI, Plots

# ## Model Setup
# The grid has length 800 km in the x-direction and 50km in the y-direction. The grid spacing is 2.5km
dx = 2500., dy = 2500., 
nx = round(Int, 800*1e3 / dx)
ny = round(Int, 50*1e3 / dy)
grid = Grid(nx=nx,ny=ny,dx=2000.,dy=2000.);

# Since we're running to steady state, we'll only do one velocity solve iteration per timestep:
solver_params = SolverParams(maxiter_picard = 1);

# and we'll take a uniform accumulation rate
params = Params( accumulation_rate = 0.3);

# ## Basal slipperiness
# We set the basal slipperiness via the model's initial conditions. In the MISMIP 3D experiment, the sliding coefficient is a Gaussian bump:
# $C = C_0 \left[1 - a  \exp (-(x- x_b)^2 / 2x_c^2  - (y - y_b)^2 / 2y_c^2) \right]

