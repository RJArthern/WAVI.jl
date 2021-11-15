```@meta
EditURL = "<unknown>/examples/melt_rate_parametrizations.jl"
```

# Melt rates

This example demonstrates the wealth of melt rate parametrizations that are included in WAVI.jl.
We produce a picture of the melt rate for the MISMIP+ steady state geometry for each of the melt rate specifications included in WAVI.jl. These are:
  * Quadratic melt rates
  * A plume melt emulator
  * PICO melt rate model
  * Binary file melt rate
You can find more info on each of these models in the Physics --> Melt Rates tab.

We also demonstrate how to add a simple melt rate model to WAVI.jl.

## Install dependencies

First let's make sure we have all required packages installed.
As well as WAVI and Plots for plotting, we're going to use the Downloads package to pull some data from a Github repository.

````@example melt_rate_parametrizations
#using Pkg
#pkg.add("https://github.com/RJArthern/WAVI.jl"), Pkg.add(Plots), Pkg.add("Downloads")
using WAVI, Plots, Downloads
````

## Setting up the  grid
First we'll make the grid, which has 2km resolution, with 320 grid points in the x-direction and 40 grid points in the y-direction:

````@example melt_rate_parametrizations
nx = 320;
ny = 40;
grid = Grid(nx=nx,ny=ny,dx=2000.,dy=2000.);
nothing #hide
````

For the bed, we'll cheat a little and use the "hard coded" form of the MISMIP+ bed (see the overdeepened bed example to see how to do this properly)

````@example melt_rate_parametrizations
bed = WAVI.mismip_plus_bed; # function definition
nothing #hide
````

## Getting the ice thickness
To save time, we're going to pull the steady state ice thickness from Github, where it is stored as a binary file.
Note that this ice thickness is the result of the "overdeepened bed" example, but with 2km resolution.

````@example melt_rate_parametrizations
fname = Downloads.download("https://github.com/alextbradley/WAVI_example_data/MISMIP_PLUS/raw/main/WAVI_ice0_2km_thick.bin");
h = Array{Float64,2}(undef, nx, ny);
read!(fname, h);
h = ntoh.(h);
nothing #hide
````

(The final line simply converts the endianness from big-endian --  the format in which this file is stored -- to little-endian.)

Now we make an InitialConditions object to store this ice thickness:

````@example melt_rate_parametrizations
initial_conditions = InitialConditions(initial_thickness = h);
nothing #hide
````

## Melt Rate Models
We'll loop over each of the melt rate models mentioned above. For each, we'll produce a map of the melt rate in the MISMIP+ geometry that we just pulled from online.

First let's make the melt rate models in turn, starting with quadratic, which is imposed on the model via a `QuadraticMeltRate` object.

````@example melt_rate_parametrizations
melt_quad = QuadraticMeltRate(γT = 0.745*1e-3)
````

The parameter $\gamma_T$ is a normalization cofficient, chosen so that the mean melt rate on the shelf is approx 10m/a. (See the WAVI Setup -->  Melt Rate Models section for information on all keyword parameters, here and for the below models.)

Next up: PICO. For this parametrization, we first have to make a mask defining where the ice front is, and then pass it when we instantiate the `PICO` object which stores the information.

````@example melt_rate_parametrizations
ice_front_mask = zeros(nx,ny);
ice_front_mask[end,:] .= 1;
melt_PICO = PICO(ice_front_mask = ice_front_mask,
                    T0 = 1.2,
                    S0 = 34.6,
                    γT = 0.87e-5,
                    nbox = 5);
nothing #hide
````

Again, $\gamma_T$ is a normalization coefficient.

Next up: plume model emulator:

````@example melt_rate_parametrizations
melt_PME = PlumeEmulator(α=1.49);
nothing #hide
````

In this case, the normalization coefficent is named $\alpha$, and has a slightly different meaning but we can use it in the same role: the control the mean melt rate on the shelf.

Finally, a binary file, in which the melt rate is read in from a file. First, we'll create such a file, which will set the melt rate to be unifrom on the shelf.

````@example melt_rate_parametrizations
isfloat = (h .< -918.0/1028.0 .* bed.(grid.xxh, grid.yyh)) #indices of floating elements
m = zeros(nx,ny);
m[isfloat] .= 10.0; #set everywhere floating to 10m/a
folder = joinpath(@__DIR__, "melt_rate_parametrizations");
isdir(folder) && rm(folder, force = true, recursive = true);
mkdir(folder) ;
out = open(joinpath(folder,"melt.bin"), "w");
write(out, m);
close(out);
nothing #hide
````

Now we can point our binary file melt rate, an instance of a `BinaryFileMeltRate` here

````@example melt_rate_parametrizations
binfile_melt = BinfileMeltRate(input_filename = joinpath(folder,"melt.bin"));
nothing #hide
````

It's useful to put these into a dictionary, so we can iterate over them:

````@example melt_rate_parametrizations
melt_rates = Dict("Quadratic" => melt_quad, "PME" => melt_PME, "PICO" => melt_PICO, "Binary file" => binfile_melt);
nothing #hide
````

## Visualizing
We'll loop over the melt rate models instantiated above. Each time, we make a model with the appropriate melt rate specifier using the 'melt_rate' keyword.
We use the `update_state!` method to bring the melt (as well as all other quantities such as grounded fraction) in line and plot the results.
We set the number of Picard iterations to 1 (we don't care about the velocity solve being completely accurate!)

````@example melt_rate_parametrizations
for (key, melt) in melt_rates
    #instantiate the model with the appropriate melt rate
    model = Model(grid = grid,
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
    #display(plt) #uncomment to show in (e.g.) VSCode
end
````

## Defining A New Melt Rate Model
In this section, we show how to specify a new melt rate model. In this example, we define a structure which can be passed to a model (via the `melt_rate` keyword) which specifies the melt rate according to the MISMIP+ experiment:
melt rate on floating cells in $0.2 \tanh((z_d - z_b)/75) \max(-100 - z_d,0)$ where $z_d$ is the ice shelf draft and $z_d - z_b$ is the cavity thickness.
There are three steps to defining a new melt rate model. First, we define a structure, which stores parameters related to the melt rate model. Note that the melt rate model does not "own" the melt rate, the `model` does (and stores it in model.fields.gh.basal_melt, see below)

````@example melt_rate_parametrizations
struct MISMIPMeltRateOne{T <: Real} <: WAVI.AbstractMeltRate
    α  :: T
    ρi :: T
    ρw :: T
end
````

In this case, a normalization coefficient $\alpha$, and the denisities of ice and ocean, $\rho_i$ and $\rho_w$, respectively.

Now we define a "constructor", a function that defines how to create one of these structures:

````@example melt_rate_parametrizations
MISMIPMeltRateOne(; α = 1.0, ρi = 918.0, ρw = 1028.0) = MISMIPMeltRateOne(α,ρi, ρw)
````

In this case, the constructor simply sets the default values for the parameters $\alpha$, $\rho_i$, and $\rho_w$.
For more complicated melt rate models, constructors might be more elaborate!

The final step is to define a function `update_melt_rate!(melt_rate::TYPE, fields, grid)` which tells WAVI how to update the melt rate in this example.
here, TYPE is the name of the structure we just made.
Note that the arguments of this function must be as mentioned here, so that the multiple dispatch capability can be leveraged! This function defines another method named `update_melt_rate!`, which sets the melt rate according to this function when the input melt rate is of type "TYPE". One of these function is defined for each of the melt rate models mentioned above.

````@example melt_rate_parametrizations
function update_melt_rate!(melt_rate::WAVI.MISMIPMeltRateOne, fields, grid)
    @unpack basal_melt, h, b  = fields.gh
    draft = -(melt_rate.ρi / melt_rate.ρw) .* h
    cavity_thickness = draft .- b
    cavity_thickness = max.(cavity_thickness, 0)
    m =  melt_rate.α .* 0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
    basal_melt[:] .= m[:]
end
````

Now we can plot the melt rate with this model:

````@example melt_rate_parametrizations
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
            title = key,
            framestyle = "box")
xlims!((420, 640))
plot!(size = (500,300))
#display(plt)
````

Finally, let's clean up the files we just made

````@example melt_rate_parametrizations
rm(joinpath(@__DIR__, "melt_rate_parametrizations"), force = true, recursive = true);
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

