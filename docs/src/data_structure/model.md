# Models

A `Model` is a WAVI.jl strcutures that contain all the information about the current state, parameters, initial conditions, and process parametrizations. An instance of a `Model` contains the following fields:
- `grid`: an instance of a WAVI.jl `Grid` object, that stores information about the numerical grid.
- `params`: an instance of a WAVI.jl `Params` object that stores phyiscal parameters that enter the model.
- `solver_params`: an instance of a WAVI.jl `SolverParams` object that stores parameters relating to the numerical solution of the governing equations (see numerical implementation).
- `initial_conditions`: an instance of a WAVI.jl `InitialConditions` object that stores initial data on the ice thickness, ice viscosity, ice temperature, and ice damage.
- `fields`: an instance of a WAVI.jl `Fields` structure that stores information relating to the current state of the model on the various grids used in the solution (see numerical implementation for more information on these grids)
- `extra_physics`: a dictionary that stores process parametrizations used in the ice sheet model (see process parametrizations for more information).

## Model Construction
A `Model` is constructed using the `Model(;<kwargs>)` constructor (a function that constructs an instance of a `Model`). Here `<kwargs>` is shorthand for keyword arguments, allow the `Model` to be configured, some are optional and others are not (see below):
```@docs
Model()
```
