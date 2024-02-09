# Simulations
The `Simulation` object is at the heart of WAVI.jl. (An instance of) a `Simulation` is constructed by passing an instance of a `Model`, an instance of a `TimesteppingParams` object, and (optionally) an instance of an `OutputParams` object, e.g.
```julia
simulation = Simulation(model = model, output_params = output_params, timestepping_params = timestepping_params)
```

```@docs
Simulation()
```
