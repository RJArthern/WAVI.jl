```@meta
EditURL = "<unknown>/examples/west_antarctica.jl"
```

# West Antarctica

This is our first real world example. We produce a map of approximate ice velocity in West Antarctica.
This example demonstrates:
  * the capability of WAVI.jl in real world examples
  * how to specify the solution mask `h_mask`, which defines which grid points are part of the ice domain.
  * how to specify a spatially variable initial ice viscosity

**NB** the files included here are intended as a full demonstration of WAVI.jl capability, but should not be interpreted as "ready to go" for scientific inquiry; please get in touch if you would like to use these files for science!

## Install dependencies

First let's make sure we have all required packages installed.
As well as WAVI and Plots for plotting, we're going to use the Downloads package to pull some data from a Github repository.

````@example west_antarctica
#using Pkg
#Pkg.add("https://github.com/RJArthern/WAVI.jl"), Pkg.add(Plots); Pkg.add("Downloads")
using WAVI, Plots, Downloads
````

## Reading in data
First, let's define our grid sizes and origin. We have a 268 x 315 with, with 3km resolution.
We're going to use 12 levels in the vertical: even though WAVI.jl is designed for the solution of depth integrated equations, it retains some information about the vertical direction e.g. in the calculation of the ice viscosity; the keyword argument `nσ`, which is passed to a `Grid` object, specifies the number of levels in the vertical. We refer to the extension of the 2D (horizontal) grid to n$\sigma$ levels in the vertical as the 3D grid, which has size nx $\times$ ny $\times$ nz.

````@example west_antarctica
nx = 268;
ny = 315;
nσ = 12;
x0 = -1792500.0;
y0 = -838500.0;
dx = 3000.0;
dy = 3000.0;
nothing #hide
````

We're going to read in the following data files:
  * h_mask     : array of zeros and ones defining the ice domain (zero corresponds to out of domain, one to in domain).
  * u_iszero   : location of grid points with zero velocity in x
  * v_iszero   : location of grid points with zero velocity in y
  * bed        : the bed elevation, which is a processed form of data from Bedmap 2.
  * h          : ice thickness , defined on the 2D grid
  * viscosity  : the initial ice viscosity, defined on the 3D grid

We need the first three of these before we can build a grid

````@example west_antarctica
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
nothing #hide
````

Before moving on to the rest, let's have a look at h_mask using the `spy` function, which indicates which entries of a matrix are non-zero:

````@example west_antarctica
#Plots.spy(hm)
````

We see a rough outline of West Antarctica! This is out ice domain. Note that indices with h_mask[i,j] = 0 are simply placeholders and do not have any physical meaning.

Next up: the bed, which will be passed to a model via the `bed_elevation` keyword article, as usual.

````@example west_antarctica
bed=Array{Float64}(undef,nx,ny);
read!(Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/WAIS/Inverse_3km_bed_clip_noNan.bin"),bed);
bed.=ntoh.(bed)
````

Let's take a look at the bed

````@example west_antarctica
plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, bed',
                    xlabel = "x (km)",
                    ylabel = "y (km)",
                    colorbar_title = "bed elevation (m)",
                    title = "West Antarctica bed elevation",
                    framestyle = "box")
````

Temperature, damage, viscosity, and thickness will be passed to a model via initial conditions (these quantities are time dependent)

````@example west_antarctica
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

initial_conditions = InitialConditions(initial_thickness = h,initial_viscosity = viscosity,initial_temperature = temp,initial_damage = damage);
nothing #hide
````

## Ice Velocity
Now we're ready to make our model, which we can then use to determine the ice velocity. All physical and solver parameters take their default values

````@example west_antarctica
model = Model(grid = grid, bed_elevation = bed,initial_conditions= initial_conditions);
nothing #hide
````

We use the `update_state!` method to bring fields (including velocity) in line with the ice thickness:

````@example west_antarctica
update_state!(model);
nothing #hide
````

Now we can visualize the ice velocity:

````@example west_antarctica
plt = Plots.heatmap(grid.xxh[:,1]/1e3, grid.yyh[1,:]/1e3, model.fields.gh.av_speed',
                    xlabel = "x (km)",
                    ylabel = "y (km)",
                    colorbar_title = "ice speed (m/yr)",
                    title = "West Antarctica ice speed",
                    framestyle = "box")
````

Ice velocities reach a maximum of approx 5km/yr on ice shelves, but are much smaller inland.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

