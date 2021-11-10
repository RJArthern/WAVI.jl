# Initial Conditions
Information on the initial state of a model is stored in an `InitialConditions` object. The following table summarizes the quantities which can be initialized in `WAVI.jl`. Each of these is set by passing an appropriately sized array to keyword argument in the constructor of (instance of an) `InitialConditions` object. Here `nx` and `ny` are the number of grid cells in the `x` and `y` directions, respectively, and `nσ` is the number of levels in the vertical. None of the keywords arguments listed below are necessary in constructing an `InitialConditions` object. (Indeed, it is not necessary to pass an `InitialConditions` object to a `Model` at all.) Quantities not specified will default to the scalar values `defualt_xxxxx` (which are stored in a `Params` object) listed below, applied everywhere.

| Keyword Argument     | Description                   | Array Size        | Default value | 
| ---------------------| ----------------------------- | ------------------- | ------------| 
| `initial_thickness`  | Ice thickness at time $t = 0$                   | `nx` x `ny`         | `default_thickness` | 
| `initial_viscosity`  |Three dimensional ice viscosity at time $t = 0$  | `nx` x `ny` x `nσ`  | `default_viscosity` | 
| `initial_temperature`|Three dimensional ice temperature at time $t = 0$  | `nx` x `ny` x `nσ`|`default_temperature` | 
| `initial_damage`     |Three dimensional ice damage at time $t = 0$  | `nx` x `ny` x `nσ`| `default_damage` | 

For example, to set the ice thickness and temperature to 500m and 265K, respectively, everywhere, we would define the following initial conditions object:
```julia
initial_conditions = InitialConditions(initial_thickness = 500.0 .* ones(nx,ny), 
                                       initial_temperature = 265.0 * ones(nx,ny,nσ))
```
We would then pass this to the constructor of a model:
```julia
model = Model(..., initial_conditions = initial_conditions)
```
In this model, the (unset) damage and viscosity would take the values `params.default_damage` and `params.default_viscosity` everywhere, respectively. 