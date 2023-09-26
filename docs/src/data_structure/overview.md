# Data Structure

WAVI.jl uses a hierachical data structure, which is shown schematically below. This page provides a brief overview of each of these structures; you can find out more information about each of these via the tabs in the sidebar.

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/build-docs/docs/src/assets/structure-flowchart.jpg" alt="" title="" width="500" height="400" /></center>
```



At the top of the hierarchy are `Simulations`. Simulations are to be ran! A `Simulation` object must be passed the sub-structures describing the following (terms in brackets are the names of the corresponding structures in `WAVI.jl`):
- Model (`Model`):  contains all the information about the current state, parameters, initial conditions, and process parametrizations (see below).
- Timestepping Parameters (`TimesteppingParams`): contains parameters relating to timestepping (e.g. timestep, number of timesteps etc)
-  Output Parameters (`OutputParams`): optional structure that contains information relativng to the outputting of solutions (what to output, when to output etc).

A `Model` structure contains sub-structures describing the following:
- Grid (`Grid`): contains information on the discretization of the model domain (e.g. number of grid cells, grid spacing etc)
- Physical parameters (`Params`): contains physical parameters (e.g. density of the ice)
- Solver parameters (`SolverParams`): contains parameters relating to the velocity solver (e.g. maximum number of iterations)
- Initial conditions (`InitialConditions`): stores initial conditions relating to the ice sheet.
- Fields (`Fields`): stores information on the current state of the model.
A `Model` also owns a `dict` named `extra_physics`, which contains information on parametrizations of physical processes used by the model. See the Parametrizations tab on the left for more information.
- Supplementary physics: this describes various other models that can be coupled to WAVI.jl. At present, this is restricted to melt rate physics, but models of damage and calving are coming soon
