# Ice Flow over a Bumpy Bed 
This is a good first two dimensional example: flow down a plane with a series of bumps superimposed. We consider how the size of the bump affects the flow speed. This example is similar in spirit to Experiment A in the Ice Sheet Model Intercomparison Exercise - Higher Order Model (ISMIP-HOM): doi: 10.5194/tc-2-95-2008
 
This example demonstrates how to
    * use WAVI.jl in two horizontal spatial dimensions with arbitrary bed shapes
    * how to interact with Grid objects.

## Install dependencies
First let's make sure we have all required packages installed.
```julia
using Pkg
Pkg.add("https://github.com/RJArthern/WAVI.jl")
Pkg.add(`Plots`)
using WAVI, Plots
```

## Basal Topography
Following ISMIP-HOM, we'll consider a bed with a series of sinusoidal oscillations with an amplitude of 500m:
$$
z_b(x,y) = -x \tan \alpha + 500 \sin (\omega x) \sin(\omega y)
$$
Here $L$ is the lengthscale of the domain, $\alpha$ is the net slope of the plane, and $\omega = 2\pi / L$ is the frequency of the bumps. We can express this bed analytically, so we write it as a function and pass it to WAVI that way
```julia
z_b(x,y; α, ω) = -x * tand(α) + 500sin(ω*x)*sin(ω*y)
```

It's useful to define our model grid here:
```julia
grid_L(; L, nx = 80,ny = 80) = Grid(nx = nx, ny = ny, dx = L/nx, dy = L/ny, y0 = 0.0, x0 = 0.0);
```
Note that `grid_L` is a function, it takes a lengthscale `L` which defines the lengthscale of the domain in both $x$ and $y$ directions (the grid is square) and optional arguments `nx` and `ny`, which set the number of grid points, with default value of 80 for both.


Let's choose a domain of 80km with 80 grid points
```julia
L = 80000.
grid80 = grid_L(L = L)
```

This grid object contains information about the location of grid points. We use this to construct an array defining the bed from the bed function `z_b` defined earlier (we don't need this step for the solver, but we do to visualize!). We'll choose the period of the bumps so that we have one peak and one trough in both directions (via the `ω` parameter)
```julia
z_b80 = z_b.(grid80.xxh,grid80.yyh; α = 0.5, ω = 2π/L );
plt = Plots.heatmap(grid80.xxh[:,1]/1e3, grid80.yyh[1,:]/1e3, z_b80, 
                        xlabel = "x (km)", 
                        ylabel = "y (km)",
                        colorbar_title = "bed depth (m)")
plot!(size = (800,800))
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//bumpy//bed.png" alt="" title="" width="600" height="600" /></center>
```

## Model Instantiation and Initial Conditions
In the ISMIP-HOM comparison, the main test is velocity along the line $y = L/4$, for various different values of $L$. Before we do that, lets look at the velocity for the example we started above with $L = 80$km. 

To begin, we create an `InitialConditions` object to prescribe the ice thickness of 1000m everywhere
```julia
initial_conditions = InitialConditions(initial_thickness = 1000. .* ones(grid80.nx, grid80.ny))
```

Now we can build our model
```julia
model80 = Model(grid = grid80, 
            bed_elevation = z_b80,
            initial_conditions = initial_conditions)
```
    
## Determining the velocity 
To bring the velocity in line with the ice thickness, we use the `update_state!` function:
```julia
update_state!(model80)
```
Now we can look at the velocity:
```julia
Plots.heatmap(model80.grid.xxh[:,1]/1e3, model80.grid.yyh[1,:]/1e3, model80.fields.gh.u', 
                        xlabel = "x (km)", 
                        ylabel = "y (km)",
                        colorbar_title = "ice velocity in x-direction (m/yr)")
plot!(size = (800,600))
```

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//bumpy//velocity_L80.png" alt="" title="" width="600" height="600" /></center>
```

## Different lengthscales
Now let's look at how the velocity along a flowline changes with the lengthscale of the domain `L`. Note that the bumps stay the same size, so as `L` increases, the aspect ratio of the bumps reduces, and we might expect they influence the flow less.

First we define the `L` values we're interested in
```julia
L_values =  [160, 80, 40, 20, 10, 5]*1.0e3;
```

We loop over these values and store the info:
```julia 
U_flowline = zeros(80, length(L_values));
grid_flowline = zeros(80, length(L_values)); # initialize storage of the the velocity (U_flowline) and x-coordinates (grid_flowline)
for (count,L) in enumerate(L_values) ;
    gridL = grid_L(L = L); #make the grid with this value of L
    z_bL = z_b.(gridL.xxh,gridL.yyh; α = 0.5, ω = 2π/L ); #make the bed with this grid
    initial_conditions = InitialConditions(initial_thickness = 1000. .* ones(80, 80)); #initial thickness of 1000m everywhere
    model = Model(grid = gridL, 
                bed_elevation = z_bL,
                initial_conditions = initial_conditions);  #build model
    update_state!(model); #get the velocity assoiciated with geometry
    grid_flowline[:, count] .= model.grid.xxh[:, round(Int, gridL.nx/4)]; #extract coordinates along line
    U_flowline[:,count] .= model.fields.gh.u[:, round(Int, gridL.nx/4)]; #get velocity along line
end
```

And make the plot:
```julia
p = plot()
#first normalize the co-ordinates
normalized_grid_flowline = zeros(size(grid_flowline))
for i = 1:(size(L_values)[1])
    plot!(grid_flowline[:,i]./L_values[i]/1e3,
    U_flowline[:,i], 
    framestyle = :box, 
    xlabel = "x/L (km)", 
    ylabel = "horizontal velocity (m/yr)",
    label = L_values[i])
end
display(p)
plot!(size = (1000,550))
```
```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/example-plots//bumpy//velocity_diffL.png" alt="" title="" width="600" height="600" /></center>
```

As expected, when the bumps have a smaller aspect ratio (smaller `L`) the flow speed is smaller.