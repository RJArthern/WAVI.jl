# Melt Rates
WAVI.jl has a number of community melt rate parametrizations/models (referred to here collectively as 'melt rates') implemented and ready to use. Details of the physics of these models can be found on the [melt rate physics page](../physics/melting.md). This page provides a description of how to endow ice sheet models with a melt rate model/parametrization, as well as the interface of each of those that are implemented in WAVI.jl.

## Adding a Melt Rate to a WAVI.jl Model
Melt rates is WAVI.jl are interfaced via `MeltRate` objects. To build an ice sheet model with a melt rate, simply build the appropriate `MeltRate` object and then pass it to the `Model` when the latter is constructed.

This is best demonstrated by example: we first build a simple ice sheet model on a 20 x 20 grid, with a bed that is 50 m below sea level:
```julia
grid = Grid(nx = 20, ny = 20)
bed_elevation = 50.0 * ones(grid.nx, grid.ny)
```
We then build our melt rate model. Let consider a simple case in which the [melt rate is assumed to depend quadratically on the thermal forcing](#Quadratic-Melt-Rate):
```julia
melt_rate = QuadraticMeltRate();
```
Finally, build the ice sheet model and pass the melt rate model
```julia
model = Model(grid = grid, bed_elevation = bed_elevation, melt_rate = melt_rate)
```


## Analytic Melt Rate Parametrizations
Documentation coming soon!

## Input File Melt Rates
Documentation coming soon!


## Quadratic Melt Rate
A `QuadraticMeltModel` -- the melt rate used to implement the [quadratic temperature melt rate parametrization](../physics/melting.md#Quadratic-Temperature-Melt-Rate-Parametrization) is constructed using the `QuadraticMeltModel(<kwargs>)` constructor. 

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

Ambient temperature and salinity profiles are passed to `PlumeEmulator` objects via the `Ta` and `Sa` keyword arguments, respectively. These must be passed as single valued functions of depth (i.e. temporal dependence in the ambient conditions is not yet supported). `Ta` and `Sa` default to the warm ambient profiles used in ISOMIP experiments (referred to as 'warm0' therein), with a lower layer of temperature 1.2∘C and salinity 34.6PSU separated from an upper layer of temperature -1∘C and salinity 33.8PSU by a pycnocline of thickness 400m, which begins at a depth of 700m below sea level.

## Plume Emulator Melt Rates
The [plume model emulator](../physics/melting.md#Plume-Emulator-Melt-Rate) of Lazeroms2018 is implemented via a `PlumeEmulator` object. Parameters used in plume model emulator melt rates are specified by keyword arguments passed to `PlumeEmulator` objects; these summarized in the following table:

| Keyword Argument   | Description                   | Units               | Default Value |
| ------------------ | ----------------------------- | ------------------- | ------------- |
| `α`                | Calibration coefficient       | Dimensionless       | 0.73          |
| `λ1`               | Liquidus slope                | ∘C                  |-0.057         |
| `λ2`               | Liquidus intercept            | ∘C                  |0.0832         |
| `λ3`               | Liquidus pressure coefficient | ∘C/m                |7.59e-4        |
| `E0`               | Entrainment coefficient       | Dimensionless       |3.6e-2         |
| `Cd`               | Drag coefficient              | Dimensionless       |2.5e-3         |
| `Γ_TS`             | Combination Stanton number    | Dimensionless       |0.0118         |
| `L`                | Latent heat of ice fusion     | J/kg                |3.35e5         |
| `c`                | Water specific heat capacity  | J/kg/∘C             |3.974e3        |
| `βs`               | Haline contraction coefficient| 1/PSU               |7.86e-4        |
| `βt`               | Thermal expansion coefficient | 1/∘C                |3.87e-5        | 
| `g`                | Gravitational acceleration    | m/s^2               |9.81           |
| `ρi`               | Ice density                   | kg/m^3              |918.0          |
| `ρw`               | Water density                 | kg/m^3              |1028.0         |

Ambient temperature and salinity profiles are passed to `PlumeEmulator` objects via the `Ta` and `Sa` keyword arguments, respectively. These must be passed as single valued functions of depth (i.e. temporal dependence in the ambient conditions is not yet supported). `Ta` and `Sa` default to the warm ambient profiles used in ISOMIP experiments (referred to as 'warm0' therein), with a lower layer of temperature 1.2∘C and salinity 34.6PSU separated from an upper layer of temperature -1∘C and salinity 33.8PSU by a pycnocline of thickness 400m, which begins at a depth of 700m below sea level.

## PICO Melt Rate Parametrization
Documentation coming soon!

## PICOp Melt Rate Parametrization
Documentation coming soon!
