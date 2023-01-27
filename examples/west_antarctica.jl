# # West Antarctica 
# 
# This is our first real world example. We produce a map of approximate ice velocity in West Antarctica.
# This example demonstrates:
#   * the capability of WAVI.jl in real world examples
#   * how to specify the solution mask `h_mask`, which defines which grid points are part of the ice domain.
#   * how to specify a spatially variable initial ice viscosity

# **NB** the files included here are intended as a full demonstration of WAVI.jl capability, but should not be interpreted as "ready to go" for scientific inquiry; please get in touch if you would like to use these files for science!

# ## Install dependencies
#
# First let's make sure we have all required packages installed. 
# As well as WAVI and Plots for plotting, we're going to use the Downloads package to pull some data from a Github repository.

#using Pkg
#Pkg.add("https://github.com/RJArthern/WAVI.jl"), Pkg.add(Plots); Pkg.add("Downloads")
using WAVI, Plots, Downloads, LinearAlgebra

BLAS.set_num_threads(1)

# ## Reading in data
# First, let's define our grid sizes and origin. We have a 268 x 315 with, with 3km resolution.
# We're going to use 12 levels in the vertical: even though WAVI.jl is designed for the solution of depth integrated equations, it retains some information about the vertical direction e.g. in the calculation of the ice viscosity; the keyword argument `nσ`, which is passed to a `Grid` object, specifies the number of levels in the vertical. We refer to the extension of the 2D (horizontal) grid to n$\sigma$ levels in the vertical as the 3D grid, which has size nx $\times$ ny $\times$ nz.
nx = 268;
ny = 315;
nσ = 12;
x0 = -1792500.0;
y0 = -838500.0;
dx = 3000.0;
dy = 3000.0;

# We're going to read in the following data files:
#   * h_mask     : array of zeros and ones defining the ice domain (zero corresponds to out of domain, one to in domain).
#   * u_iszero   : location of grid points with zero velocity in x
#   * v_iszero   : location of grid points with zero velocity in y
#   * bed        : the bed elevation, which is a processed form of data from Bedmap 2.
#   * h          : ice thickness , defined on the 2D grid
#   * viscosity  : the initial ice viscosity, defined on the 3D grid 

# We need the first three of these before we can build a grid
h_mask=Array{Float64}(undef,nx,ny);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_h_mask_clip.bin"),h_mask);
hm = ntoh.(h_mask);
hm = map.(Bool, round.(Int, hm));

u_iszero=Array{Float64}(undef,nx+1,ny);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_uiszero_clip.bin"),u_iszero);
u_iszero.=ntoh.(u_iszero);
u_iszero = map.(Bool, round.(Int, u_iszero));

v_iszero=Array{Float64}(undef,nx,ny+1);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_viszero_clip.bin"),v_iszero);
v_iszero.=ntoh.(v_iszero);
v_iszero = map.(Bool, round.(Int, v_iszero));

grid = Grid(nx = nx, ny = ny, nσ = nσ, x0 = x0, y0 = y0, dx = dx, dy = dy, h_mask = hm, u_iszero = u_iszero, v_iszero = v_iszero);

# Before moving on to the rest, let's have a look at h_mask using the `spy` function, which indicates which entries of a matrix are non-zero:
#Plots.spy(hm)
# We see a rough outline of West Antarctica! This is out ice domain. Note that indices with h_mask[i,j] = 0 are simply placeholders and do not have any physical meaning.

# Next up: the bed, which will be passed to a model via the `bed_elevation` keyword article, as usual.
bed=Array{Float64}(undef,nx,ny);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_bed_clip_noNan.bin"),bed);
bed.=ntoh.(bed)

# Let's take a look at the bed
plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, bed', 
                    xlabel = "x (km)", 
                    ylabel = "y (km)",
                    colorbar_title = "bed elevation (m)",
                    title = "West Antarctica bed elevation",
                    framestyle = "box")


# Temperature, damage, viscosity, and thickness will be passed to a model via initial conditions (these quantities are time dependent)
temp=Array{Float64}(undef,nx,ny,nσ);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_3Dtemp_clip_noNan.bin"),temp);
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_damage3D_clip_noNan.bin"),damage);
damage.=ntoh.(damage)

h=Array{Float64}(undef,nx,ny);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_thickness_clip_noNan.bin"),h);
h.=ntoh.(h)

viscosity=Array{Float64}(undef,nx,ny,nσ);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_viscosity3D_clip_noNan.bin"),viscosity);
viscosity.=ntoh.(viscosity);

solver_params=SolverParams(maxiter_picard=20)

initial_conditions = InitialConditions(initial_thickness = h,initial_viscosity = viscosity,initial_temperature = temp,initial_damage = damage);

parallel_spec = BasicParallelSpec()

# ## Ice Velocity
# Now we're ready to make our model, which we can then use to determine the ice velocity. All physical and solver parameters take their default values
serial_model = Model(grid = grid, bed_elevation = bed,initial_conditions= initial_conditions, parallel_spec = parallel_spec, solver_params=solver_params);

# We use the `update_state!` method to bring fields (including velocity) in line with the ice thickness:
update_state!(serial_model)

# Now we can visualize the ice velocity:plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, log10.(serial_model.fields.gh.av_speed'), 
plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, log10.(serial_model.fields.gh.av_speed'), 
                    xlabel = "x (km)", 
                    ylabel = "y (km)",
                    colorbar_title = "ice speed (m/yr)",
                    title = "West Antarctica ice speed",
             #       clims=(-0.02,0.02),
             #       xlims=(-1400,-1300),
             #       ylims=(-500,-400),
                    framestyle = "box")


# Now try the same velocity solve with a shared memory parallel implementation                     
saved_serial_model=deepcopy(serial_model);

parallel_model = deepcopy(serial_model)
parallel_model = @set parallel_model.solver_params = SolverParams(maxiter_picard=5)
parallel_model = @set parallel_model.parallel_spec = SharedMemorySpec(ngridsx=2,ngridsy=3,overlap=1,damping=0.0,niterations=1)

cols=cgrad(:roma,rev=true);
nSchwarzIterations=1
for iSchwarzIteration = 1:nSchwarzIterations
       update_state!(parallel_model)


       # Now we can visualize the ice velocity:
       plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, log10.(parallel_model.fields.gh.av_speed'), 
                     c = cols, 
                     xlabel = "x (km)", 
                     ylabel = "y (km)",
                     colorbar_title = "ice speed (m/yr)",
                     title = "West Antarctica ice speed",
              #       clims=(-0.2,0.2),
              #       xlims=(-1600,-1500),
              #       ylims=(-700,-600),
                     framestyle = "box")
       display(plt)

       plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, log10.(parallel_model.fields.gh.ηav'), 
       c = cols,
       xlabel = "x (km)", 
       ylabel = "y (km)",
       colorbar_title = "ice speed (m/yr)",
       title = "West Antarctica ice speed",
#       clims=(-0.2,0.2),
#       xlims=(-1600,-1500),
#       ylims=(-700,-600),
       framestyle = "box")
       display(plt)

       plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, log10.(parallel_model.fields.gh.av_speed'./serial_model.fields.gh.av_speed'), 
       c = cols,
       xlabel = "x (km)", 
       ylabel = "y (km)",
       colorbar_title = "ice speed (m/yr)",
       title = "West Antarctica ice speed",
#       clims=(-0.2,0.2),
#       xlims=(-1600,-1500),
#       ylims=(-700,-600),
       framestyle = "box")
       display(plt)

# Ice velocities reach a maximum of approx 5km/yr on ice shelves, but are much smaller inland.

       plt = Plots.heatmap(grid.xxu[:,1]/1e3, grid.yyu[1,:]/1e3, log10.(abs.(parallel_model.fields.gu.residual')), 
                     c = cols,
                     xlabel = "x (km)", 
                     ylabel = "y (km)",
                     colorbar_title = "ice speed (m/yr)",
                     title = "West Antarctica ice speed",
              #       clims=(-0.02,0.02),
              #       xlims=(-1400,-1300),
              #       ylims=(-500,-400),
                     framestyle = "box")
       display(plt)

       plt = Plots.heatmap(grid.xxv[:,1]/1e3, grid.yyv[1,:]/1e3, log10.(abs.(parallel_model.fields.gv.residual')), 
                     c=cols,
                     xlabel = "x (km)", 
                     ylabel = "y (km)",
                     colorbar_title = "ice speed (m/yr)",
                     title = "West Antarctica ice speed",
              #       clims=(-0.02,0.02),
              #       xlims=(-1400,-1300),
              #       ylims=(-500,-400),
                     framestyle = "box")
       display(plt)

 end

