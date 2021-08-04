# Melt Rate Models
WAVI.jl has a number of community melt rate parametrizations/models implemented and ready to use. Details of these models can be found on the [melt rate physics page](../physics/melting.md). This page provides a description of how to endow ice sheet models with a melt rate model/parametrization, as well as the interface of each of those that are implemented in WAVI.jl.

## Adding a Melt Rate Parametrization to a WAVI.jl Model
To add a melt rate parametrization to a WAVI.jl model, simply build your (ice sheet) `Model` and particular `MeltRateModel`, and use the `add_melt_rate_model!` function to glue them together. 

This is best demonstrated by example: we first build a simple ice sheet model on a 20 x 20 grid, with a bed that is 50 m below sea level:
```julia
grid = Grid(nx = 20, ny = 20)
bed_elevation = 50.0 * ones(grid.nx, grid.ny)
model = Model(grid = grid, bed_elevation = bed_elevation)
```
We then build our melt rate model. In this example, a simple [analytic parametrization of the melt rate](#Analytic-Parametrization) that specifies the melt rate according to some user defined function (in this case the maximum of zero and the ice draft minus 100):
```julia
function melt_function())
    draft = -(918.0 / 1028.0) * h
    m =  max.((-100 .- draft), 0)
    return m
end
arguments = (h = model.fields.gh.h,)
melt_rate_model = AnalyticMeltRate(melt_rate_function = melt_function, 
                                function_arguments = arguments)
```
Finally, the ice sheet model is endowed with the melt rate mode:
```julia
add_melt_rate_model!(model, melt_rate_model)
```
Simulations that use this model will have evolve with melt rates specified according to the function `melt_function`.

## Analytic Melt Rate Parametrizations

## Input File Melt Rates

## Quadratic Temperature Melt Rate Model
A `QuadraticMeltModel` -- the melt model used to implement the [quadratic temperature melt rate parametrization](../physics/melting.md#Quadratic-Temperature-Melt-Rate-Parametrization) is constructed using the `QuadraticMeltModel(<kwargs>)` constructor. 

A `QuadraticMeltModel` accepts the following keyword arguments:
- `h` (required): array of ice thickness values at grid points
- `melt_partial_cell` (default: `false`): specify whether to apply melt to partially floating cells or not.
- `lambda_1` (default: `-0.057`): liquidus slope
- `lambda_2` (default: `0.0832`): liquidus intercept
- `lambda_3` (default: `7.59e-4`): liquidus pressure coefficient.
- `gamma_T` (default `99.32e-5`): heat exchange velocity. Note that `gamma_T` is typically used as a tuning parameter, the default value is the tuned value from [Favier2019](@cite).
- `L` (default `3.35e-5`): latent heat of fusion
- `rho_s` (default `1028.0`): sea water density
- `rho_i` (default `918.0`): ice density
- `c_p` (default 3974.0): specific heat capacity of ocean
- `S_0`: far-field practical salinity profile (units: psu), passed to the constructor as a function of depth. The default is the 'warm0' profile used in [Favier2019](@cite).
- `T_0`: far-field potential temperature profile (units: ${}^\circ$C), passed to the constructor as a function of depth. The default is the 'warm0' profile used in [Favier2019](@cite).



## Plume Emulator Melt Rate Parametrization

## PICO Melt Rate Parametrization

## PICOp Melt Rate Parametrization