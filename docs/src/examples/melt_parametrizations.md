## Melt rates parametrizations

This example demonstrates the how to use melt rate parametrizations in WAVI.jl. We first demonstrate the parametrizations included in WAVI.jl, by producing a map of the melt rate for each of these in the MISMIP+ steady state geometry. (MISMIP+ is the latest ice sheet model intercomparison exercise, for more info see doi:10.5194/tc-14-2283-2020)
Explicitly, these parametrizations are:
   * Quadratic melt rates
   * A plume melt emulator
   * PICO melt rate model
   * Binary file melt rate
You can find more info on each of these models in the Physics --> Melt Rates tab.
 
Secondly, we demonstrate how to add a simple melt rate model to WAVI.jl.


## Install dependencies

First let's make sure we have all required packages installed. As well as WAVI and Plots for plotting, we're going to use the Downloads package to pull some data from a Github repository.
```julia
using Pkg
pkg.add("https://github.com/RJArthern/WAVI.jl"), 
Pkg.add("Plots"), Pkg.add("Downloads")
using WAVI, Plots, Downloads
```

## Model Setup
First we'll make the grid. The MISMIP+ grid has length 640km and 80km in the x and y directions, respectively. We're going to choose 2km resolution in both directions, i.e. 320 grid points in the x-direction and 40 grid points in the y-direction.
```julia
nx = 320;
ny = 40;
grid = Grid(nx=nx,ny=ny,dx=2000.,dy=2000.);
```

For the bed, we'll cheat a little and use the "hard coded" form of the MISMIP+ bed (see the "MISMIP+" example to see how to define this bed properly)

```julia
bed = WAVI.mismip_plus_bed # function definition
```

To allow us to focus on the melt rate parametrizations, we're not going to run a simulation to define the steady state ice thickness, but rather pull it from Github, where it is stored as a binary file. (Note that this ice thickness is the result of the "MISMIP+" example, but with 2km resolution.)
```julia
fname = Downloads.download("https://github.com/alextbradley/WAVI_example_data/raw/main/MISMIP_PLUS/WAVI_ice0_2km_thick.bin");
h = Array{Float64,2}(undef, nx, ny);
read!(fname, h);
h = ntoh.(h);
# (The final line simply converts the endianness from big-endian --  the format in which this file is stored -- to little-endian.)
```

Now we make an InitialConditions object to store this ice thickness
```julia
global initial_conditions = InitialConditions(initial_thickness = h);
```
(We make this a global variable because we'll use it in a loop later.)

## Melt Rate Models  
We'll loop over each of the melt rate models mentioned above. For each, we'll produce a map of the melt rate in the MISMIP+ geometry that we just downloaded.

First let's make the melt rate models in turn, starting with quadratic, which is imposed on the model via a `QuadraticMeltRate` object
```julia
melt_quad = QuadraticMeltRate(γT = 0.745*1e-3)
```
The parameter $\gamma_T$ is a normalization cofficient. The value here is chosen so that the mean melt rate on the shelf is approx 10m/a. (See the "WAVI Setup" -->  "Melt Rate Models" section for information on all keyword parameters, here and for the below models.)

Next up: PICO. For this parametrization, we first have to make a mask defining where the ice front is, and then pass it as a keyword argument when we instantiate the `PICO` melt rate object.
```julia
ice_front_mask = zeros(nx,ny);
ice_front_mask[end,:] .= 1; #this tells julia that the ice front is at the downstream end of the domain
melt_PICO = PICO(ice_front_mask = ice_front_mask, 
                    T0 = 1.2, 
                    S0 = 34.6, 
                    γT = 0.87e-5, 
                    nbox = 5);
```
Again, $\gamma_T$ is a normalization coefficient.

Next up: plume model emulator:
```julia
melt_PME = PlumeEmulator(α=1.49);
```
In this case, the normalization coefficent is named $\alpha$, and has a slightly different meaning to $\gamma_T$ above, but we use it in the same role: to the set the mean melt rate on the shelf.

Finally, a binary file melt rate, in which the melt rate is read in from a binary file. First, we'll create such a file, which will set the melt rate to be uniform on the shelf.
```julia
isfloat = (h .< -918.0/1028.0 .* bed.(grid.xxh, grid.yyh)) #indices of floating elements
m = zeros(nx,ny);
m[isfloat] .= 10.0; #set everywhere floating to 10m/a
folder = joinpath(@__DIR__, "melt_rate_parametrizations");
isdir(folder) && rm(folder, force = true, recursive = true);
mkdir(folder) ;
out = open(joinpath(folder,"melt.bin"), "w");
write(out, m);
close(out);
```

Now we instantiate out melt rate opject, pointing it to this file we just created:
```julia
binfile_melt = BinfileMeltRate(input_filename = joinpath(folder,"melt.bin"))
```

It's useful to put these into a dictionary, so we can iterate over them. Again, we use a global variable to facilitate this
```julia
global melt_rates = Dict("Quadratic" => melt_quad, "PME" => melt_PME, "PICO" => melt_PICO, "Binary file" => binfile_melt)
```

## Visualization
We'll loop over the melt rate models instantiated above. Each time, we make a WAVI.jl `Model` with the appropriate melt rate parametrization specified via the `melt_rate` keyword.  We use the `update_state!` method to bring the melt (as well as all other quantities, such as the grounded fraction) in line with the specified thickness, and then plot the melt rate. (We also set the number of Picard iterations to 1 -- each time we do an `update_state!` we perform a velocity solve; in this example, we don't care about the velocity solve being completely accurate, so we'll do the minimum number of iterations in this process!)
```julia
for (key, melt) in melt_rates
    #instantiate the model with the appropriate melt rate
    local model = Model(grid = grid,
            bed_elevation = bed, 
            initial_conditions = initial_conditions,
            solver_params = SolverParams(maxiter_picard=1),
            melt_rate = melt);

    #bring quantities in line with the thickness
    update_state!(model);

    #extract the melt rate, remove any grounded entries and saturate the melt rate to 50 m/a
    mcopy= deepcopy(model.fields.gh.basal_melt)
    mcopy[model.fields.gh.grounded_fraction .== 1.] .= NaN
    msat = deepcopy(mcopy)
    msat[msat .> 50] .= 50

    #plot the melt rate
    plt = Plots.heatmap(model.grid.xxh[:,1]/1e3, model.grid.yyh[1,:]/1e3, msat', 
                        xlabel = "x (km)", 
                        ylabel = "y (km)",
                        colorbar_title = "melt rate (m/yr)",
                        title = key,
                        framestyle = "box")
    xlims!((420, 640))
    plot!(size = (500,300))
    display(plt) #uncomment to show in (e.g.) VSCode
end
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//melt_parametrizations//quadratic.png" alt="" title="" width="600" height="600" /></center>
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//melt_parametrizations//pico.png" alt="" title="" width="600" height="600" /></center>
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//melt_parametrizations//plume.png" alt="" title="" width="600" height="600" /></center>
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//melt_parametrizations//binary.png" alt="" title="" width="600" height="600" /></center>
```

As a sanity check, the binary file melt rate has the same (10m/a)melt rate over the whole shelf. The PICO parametrization, which divides the shelf up into discrete chunks, has a corresponding banded structure, with highest melt rates at the grounding line (note the different colourbar limits on the various plots!). The quadratic melt rate parametrizations similarly has the highest melt rate near the grounding line, but drops off with distance from the grounding line much quicker than the PICO parametrization. These plots can be compared to corresponding results for the NEMO ocean model  (Favier et al. 2019 doi:10.5194/gmd-12-2255-2019)

## Defining A New Melt Rate Model
 In this section, we show how to specify a new melt rate model. NB: This is experimental and still in development. There are four steps:
# * Create a file to store code
# * Define the appropriate structure, which stores information required to prescribe the melt rate
# * Define a constructor of this structure
# * Write a function to update the melt rate appropriately, and export this
# This is quite abstract, so let's do an example. We'll create a melt rate model which sets the melt rate as it is specified in the MISMIP+ experiment: 
# the melt rate on floating cells is $0.2 \tanh((z_d - z_b)/75) \max(-100 - z_d,0)$, where $z_d$ is the ice shelf draft and $z_d - z_b$ is the cavity thickness. 
# We'll follow the steps above: first, we create a file to store the code. For this example, we've already create the file, you can see it at it at "src/MeltRate/mismip_melt_rate.jl".

# Next we define a structure, which stores parameters related to the melt rate model. Note that the melt rate model does not "own" the melt rate, the `model` does (and stores it in model.fields.gh.basal_melt, see below)
""" 
struct MISMIPMeltRateOne{T <: Real} <: AbstractMeltRate 
    α  :: T
    ρi :: T
    ρw :: T
end
"""
# In this case, a normalization coefficient $\alpha$, and the denisities of ice and ocean, $\rho_i$ and $\rho_w$, respectively. 

# Now we define our "constructor", a function that defines how to create one of these structures:
""" 
MISMIPMeltRateOne(; α = 1.0, ρi = 918.0, ρw = 1028.0) = MISMIPMeltRateOne(α,ρi, ρw)
""" 
# In this case, the constructor simply sets the default values for the parameters $\alpha$, $\rho_i$, and $\rho_w$. 
# For more complicated melt rate models, constructors might be more elaborate!

# The final step is to define a function `update_melt_rate!(melt_rate::TYPE, fields, grid)` which tells WAVI how to update the melt rate in this example.
# here, TYPE is the name of the structure we just made. 
# Note that the arguments of this function must be as mentioned here, so that the multiple dispatch capability can be leveraged! This procedure defines another method named `update_melt_rate!`, which sets the melt rate according to this function when the input melt rate is of type "TYPE". One of these function is defined for each of the melt rate models mentioned above.
""" 
function update_melt_rate!(melt_rate::MISMIPMeltRateOne, fields, grid) 
    draft = -(melt_rate.ρi / melt_rate.ρw) .* fields.gh.h
    cavity_thickness = draft .- fields.gh.b
    cavity_thickness = max.(cavity_thickness, 0)
    m =  melt_rate.α .* 0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
    fields.gh.basal_melt[:] .= m[:]
end
""" 

# Finally, we tell point WAVI to this code by adding `include("mismip_melt_rate.jl")` to the file "src/MeltRate/MeltRate.jl", and add this structure to the export section in "src/WAVI.jl" file. 

# Now we can create a model which takes this melt rate andplot the result:
model = Model(grid = grid,
            bed_elevation = bed, 
            initial_conditions = initial_conditions,
            solver_params = SolverParams(maxiter_picard=1),
            melt_rate = MISMIPMeltRateOne());

update_state!(model);
m = deepcopy(model.fields.gh.basal_melt)
m[model.fields.gh.grounded_fraction .== 1.] .= NaN
msat = deepcopy(m)
msat[msat .> 50] .= 50

plt = Plots.heatmap(model.grid.xxh[:,1]/1e3, model.grid.yyh[1,:]/1e3, msat', 
            xlabel = "x (km)", 
            ylabel = "y (km)",
            colorbar_title = "melt rate (m/yr)",
            title = "MISMIP melt rate",
            framestyle = "box")
xlims!((420, 640))
plot!(size = (500,300))
#display(plt)

# Hopefully this example demonstrates clearly the procedure for adding melt rate models to WAVI.jl. If there are any questions, don't hesistate to get in touch (see the "Contact Us" tab)

# Finally, let's clean up the files we just made
rm(joinpath(@__DIR__, "melt_rate_parametrizations"), force = true, recursive = true);
