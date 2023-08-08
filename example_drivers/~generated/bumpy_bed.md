```@meta
EditURL = "<unknown>/examples/bumpy_bed.jl"
```

# Ice Flow over a Bumpy Bed

This is the first two dimensional example: flow down a plane with bumps with different length scales.
This example is similar in spirit to Experiment A in the Ice Sheet Model Intercomparison Exercise - Higher Order Model (ISMIP-HOM): doi: 10.5194/tc-2-95-2008
This example demonstrates: How to use WAVI.jl in plan view (two horizontal spatial dimensions) with arbitrary bed shapes, and how to interact with Grid objects.

## Install dependencies

First let's make sure we have all required packages installed.

using Pkg
Pkg.add("https://github.com/RJArthern/WAVI.jl"), Pkg.add(Plots)

````@example bumpy_bed
using WAVI, Plots
````

## Basal Topography

Following ISMIP-HOM, let's take a bed with a series of sinusoidal oscillations with an amplitude of 500m:
$$
z_b(x,y) = -x \tan \alpha + 500 \sin (\omega x) \sin(\omega y)
$$
Here $L$ is the lengthscale of the domain, $\alpha$ is the net slope of the plane, and $\omega = 2\pi / L$ is the frequency of the bumps.
Let's define this bed as a function and take a look at it:

````@example bumpy_bed
z_b(x,y; α, ω) = -x * tand(α) + 500sin(ω*x)*sin(ω*y);
nothing #hide
````

It's useful to define our model grid here. We can do this for an arbitary domain length L.

````@example bumpy_bed
grid(; L, nx = 80,ny = 80) = Grid(nx = nx, ny = ny, dx = L/nx, dy = L/ny, y0 = 0.0, x0 = 0.0);
nothing #hide
````

The final two arguments simply set the origin of the co-ordinate system. This grid has 80 grid points in each direction by default

Let's choose a domain of 80km with 80 grid points

````@example bumpy_bed
L = 80000.;
grid80 = grid(L = L);
nothing #hide
````

This grid object contains information about the location of grid points. We use this to construct and array defining the bed:

````@example bumpy_bed
z_b80 = z_b.(grid80.xxh,grid80.yyh; α = 0.5, ω = 2π/L );
plt = Plots.heatmap(grid80.xxh[:,1]/1e3, grid80.yyh[1,:]/1e3, z_b80,
                        xlabel = "x (km)",
                        ylabel = "y (km)",
                        colorbar_title = "bed depth (m)")
plot!(size = (800,400))
#display(plt)
````

## Model Instantiation and Initial Conditions
In the ISMIP-HOM comparison, the main test is velocity along the line $y = L/4$, for various different values of $L$.
Before we do that, lets look at the velocity for the example we started above with $L = 80$km.

To begin, we create an `InitialConditions` object to prescribe the ice thickness of 1000m everywhere

````@example bumpy_bed
initial_conditions = InitialConditions(initial_thickness = 1000. .* ones(grid80.nx, grid80.ny));
nothing #hide
````

Now we can build our model

````@example bumpy_bed
model80 = Model(grid = grid80,
            bed_elevation = z_b80,
            initial_conditions = initial_conditions);
nothing #hide
````

## Determining the velocity
To bring the velocity in line with the ice thickness, we have to use the `update_state!` function

````@example bumpy_bed
update_state!(model80);
nothing #hide
````

Let's have a look at the velocity component in the x-direction

````@example bumpy_bed
Plots.heatmap(model80.grid.xxh[:,1]/1e3, model80.grid.yyh[1,:]/1e3, model80.fields.gh.u',
                        xlabel = "x (km)",
                        ylabel = "y (km)",
                        colorbar_title = "ice velocity in x-direction (m/yr)")
plot!(size = (1000,550))
````

## Different Lengthscales
Now let's look at how the velocity along a flowline changes with L.

````@example bumpy_bed
L_values =  [160, 80, 40, 20, 10, 5]*1.0e3;
nothing #hide
````

We loop over these values and store the info:

````@example bumpy_bed
U_flowline = zeros(80, length(L_values));
grid_flowline = zeros(80, length(L_values));
for (count,L) in enumerate(L_values) ;
    gridL = grid(L = L);
    z_bL = z_b.(gridL.xxh,gridL.yyh; α = 0.5, ω = 2π/L );
    initial_conditions = InitialConditions(initial_thickness = 1000. .* ones(80, 80));
    model = Model(grid = gridL,
                bed_elevation = z_bL,
                initial_conditions = initial_conditions);
    update_state!(model);
    grid_flowline[:, count] .= model.grid.xxh[:, round(Int, gridL.nx/4)];
    U_flowline[:,count] .= model.fields.gh.u[:, round(Int, gridL.nx/4)];
end

#And make the plot in a single command
plot(grid_flowline/1e3,
    U_flowline,
    layout = (2,3),
    framestyle = :box,
    xlabel = "x (km)",
    ylabel = "horizontal velocity (m/yr)",
    label = :none,
    title =  ["160km" "80km" "40km" "20km" "10km" "5km"])
plot!(size = (1000,550))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

