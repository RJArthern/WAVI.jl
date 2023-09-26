<!-- Title -->
<h1 align="center">
  ☃️🏔️❄️ WAVI.jl ❄️🏔️☃️
</h1>

WAVI (Wavelet-based Adaptive-grid Vertically-integrated Ice-model) is a fast and friendly ice sheet model, written in Julia. WAVI documentation can be found [here](https://rjarthern.github.io/WAVI.jl/).

## Contents

* [Installation instructions](#installation-instructions)
* [Running your first model](#running-your-first-model)
* [Getting help](#getting-help)
* [Contributing](#contributing)
* [Credits](#credits)

## Installation instructions
You can install the latest version of WAVI using Julia's in-build package manager:
```julia
julia>using Pkg
julia>Pkg.add(PackageSpec(url="https://github.com/RJArthern/WAVI.jl.git", rev = "main"))
```
Note that WAVI requires Julia v1.5 or newer.

Updating WAVI is also achieved using the package manager
```julia
julia>using Pkg
julia>Pkg.update("WAVI"))
```
Note that updating should be done with care as WAVI is still developing rapidly; while we aim to keep breaking changes to a minimum, this cannot be guaranteed at present.

## Running your first model
Let's run the MISMIP+ experiment (http://www.climate-cryosphere.org/activities/targeted/153-misomip/1412-mismip-plus), the latest ice sheet model intercomparison experiment. We'll use a grid with 80x10 grid points and 8km resolution in both dimensions, and 4 vertical levels. We'll run the model to steady state for 10000 years with a timestep of 0.5 years. Since we're only interested in the steady state, we'll speed up the code by only doing one iteration of the velocity solve per timestep:
```julia
using WAVI 
grid = Grid(nx = 80, ny = 10, nσ = 4, dx = 8000., dy = 8000., u_iszero = ["north"], v_iszero = ["east", "west"])
bed = WAVI.mismip_plus_bed 
params = Params(accumulation_rate = 0.3)
solver_params =  SolverParams(maxiter_picard = 1)
model = Model(grid = grid,bed_elevation = bed, params = params)
timestepping_params = TimesteppingParams(dt = 0.5, end_time = 10000.)
simulation = Simulation(model = model, timestepping_params = timestepping_params)
run_simulation!(simulation)
```
It's as easy as that: entry into the state of the art ice sheet model intercomparison in nine lines of code 😎

## Getting help

## Contributing

## Credits
This package was initiated by Rob Arthern (https://github.com/RJArthern) and is currently maintained by Rob and Alex Bradley (https://github.com/alextbradley)
