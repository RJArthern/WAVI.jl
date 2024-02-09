# Timestepping Parameters

## Overview
A `TimesteppingParams` structure is used to pass information relating to [timestepping](../numerical_procedure/numerical_procedure.md) to the simulation. The following parameters as specified by passing appropriate keyword arguments the `TimesteppingParams` constructor:
- `niter0`: the iteration number at which the simulation. Set `niter0` to `0` to initialize a 'clean' simulation, or `niter0` to a positive integer to (attempt to) pickup from the a corresponding checkpoint (see [Checkpoints and Pickups](#Checkpoints-and-Pickups) below)
- `dt`: the simulation timestep
- `end_time`: the clock time at which the simulation should terminate
- `n_iter_total`: the total number of timesteps to be performed. **NB**: you must specify at least one of `end_time` and `n_iter_total` (the simulation must know when it is going to finish!); specifying both is possible, but they must be compatible (i.e the `end_time` must equal `n_iter_total * dt`).
- `chkpt_freq` and `pchkpt_freq`: the model time interval between outputting temporary and permanent checkpoints, respectively (see [Checkpoints and Pickups](#Checkpoints-and-Pickups))
- `step_thickness`: a flag which, when turned off (`stepthickness = false`) turns of thickness updates when timestepping. This is an experimental feature, necessary for coupling WAVI to the MITgcm (see the [MITgcm coupling](../mitgcm_coupling.md) for more info)

## Constructor
An instance of a `TimesteppingParams` is constucted using the `TimesteppingParams(;<kwargs>)` constructor:
```@docs
TimesteppingParams()
```

## Checkpoints and Pickups
Large simulations are computationally expensive, and may take a long time to run. To permit simulations to run for longer than maximum runtime limits which are imposed on many machines, WAVI.jl is equipped with a checkpoint-pickup system that allows the state to be outputted frequently, and the simulation to be picked again from that point. 

### Checkpoints
Checkpoints contain a snapshot of *all* of the information on the current state of the simulation. As a result, they may be large and thus outputting checkpoints very frequently is discouraged, except for when debugging simulations which are ran for only a few timesteps (see [Simulation Tips](../simulation_tips.md))

Both temporary and permanent checkpoints are output in `jld2` format. The path of both of these is specified by setting the `chkpt_path` argument in the TimesteppingParams [constructor](#Constructor).

Temporary checkpoints are outputted every `chkpt_freq` years (model time). These checkpoints are rolling: first checkpoint 'A' (named `chkptA.jld`) is written at time `chkpt_freq`, then checkpoint 'B' (named `chkptB.jld`) is output at a time `2 * chkpt_freq`. Checkpoint 'A' is overwritten at time `3 * chkpt_freq`, and so on.

Similarly, permanent checkpoints are outputted every `pchkpt_freq` model years.  These are output sequentially (i.e. not overwritten) as `pchkpt_xxxxxxxxx.jld2` where `xxxxxxxxxx` is the iteration number, prepended by zeros to be ten digits long (e.g. the filename of a checkpoint outputted at timestep number one thousand is `chkpt_0000001000.jld`).

Note that the number of timesteps, rather than the model time, is used to determine when to output the solution: permanent checkpoints are output every `floor(pchkpt_freq / dt)` timesteps (and similarly for `chkpt_freq`), so the actual output frequency may not exaclty match that specified by the parameter. 

### Pickups
Continuining from the state stored in a checkpoint is initiatied by setting the `niter0` parameter to correspond to the number of the timestep to be picked up from. For example, to continue the simulation from a checkpoint outputted at timestep number one thousand, set `niter0 = 1000` in the `TimesteppingParams` constructor.  

WAVI.jl only currently supports pickups from permanent checkpoints. The permanent checkpoint from which a pickup is initiated must be in the local directory.

!!! note
    **After a pickup, WAVI.jl uses the model that was specified in the original simulation that produced the checkpoint**. However, to retain the simulation structure, a `model` must still be passed to `simulation` Similarly, **after picking up, `output_params` will be as specified in the original simulation that produced the checkpoint**. To force a simulation to update to the model and/or output parameters to those just passed, use the `pickup_model_update_flag` and/or `pickup_output_update_flag` flags in the simulation constructor (see [Simulations](simulations.md)). **`TimesteppingParams` are always updated after a pickup**.  

