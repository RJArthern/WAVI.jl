# Timestepping Parameters

## Overview
A `TimesteppingParams` structure is used to pass information relating to [timestepping](../numerical_procedure/numerical_procedure.md) to the simulation. The following parameters are specified by passing appropriate keyword arguments to the `TimesteppingParams` constructor:

- `niter0`: the iteration number at which the simulation starts. Set `niter0` to `0` to initialize a new simulation, or to a positive integer to pick up from a checkpoint written at that iteration (see [Checkpoints and pickups](#Checkpoints-and-pickups) below).
- `dt`: the simulation timestep
- `end_time`: the clock time at which the simulation should terminate
- `n_iter_total`: the total number of timesteps to be performed. **NB**: you must specify at least one of `end_time` and `n_iter_total` (the simulation must know when it is going to finish!); specifying both is possible, but they must be compatible (i.e. `end_time` must equal `n_iter_total * dt`).
- `chkpt_freq`: the model time interval between writing permanent checkpoints (see [Checkpoints and pickups](#Checkpoints-and-pickups)). Set to `Inf` (the default) to disable checkpointing.
- `chkpt_path`: directory in which checkpoint files are stored when `chkpt_path` is set explicitly (see [Checkpoint location](#Checkpoint-location) below).
- `step_thickness`: a flag which, when turned off (`step_thickness = false`), turns off thickness updates when timestepping. This is an experimental feature, necessary for coupling WAVI to the MITgcm (see the [MITgcm coupling](../mitgcm_coupling.md) for more info).

The deprecated keyword `pchkpt_freq` is no longer accepted; use `chkpt_freq` instead.

## Constructor
An instance of `TimesteppingParams` is constructed using the `TimesteppingParams(; kwargs...)` constructor:

```@docs
TimesteppingParams()
```

## Checkpoints and Pickups
Large simulations are computationally expensive, and may take a long time to run. To permit simulations to run for longer than maximum runtime limits which are imposed on many machines, WAVI.jl is equipped with a checkpoint-pickup system that allows the state to be outputted frequently, and the simulation to be picked again from that point.

### Checkpoints
Checkpoints contain a snapshot of the evolving ice state and the simulation clock. They can be large, so writing them very frequently is discouraged except for short test runs (see [Simulation tips](../simulation_tips.md)).

Permanent checkpoints are written in `jld2` format. Each file is named `Chkpt_NNNNNNNNNN.jld2`, where `NNNNNNNNNN` is the iteration number zero-padded to ten digits (for example, iteration `1000` gives `Chkpt_0000001000.jld2`).

The snapshot stores:

- **Clock:** `n_iter`, `time`, `ref_time`
- **Grid geometry:** `nx`, `ny`, `nσ`, `dx`, `dy`, `x0`, `y0`
- **Bed**
- **`InitialConditions` arrays:** thickness, grounded fraction, velocities, viscosity, temperature, damage, strain history, basal water, hydraulic potential, effective pressure, depth-averaged temperature, and `preBfactor`
- **Basal drag** `β` (not part of `InitialConditions`, but needed for restart)

Wavelets and other physics objects are not in the file.

Checkpoints are written when all of the following hold:

- `chkpt_freq` is not `Inf`
- the current iteration index is greater than zero
- `mod(n_iter, n_iter_chkpt) == 0`, where `n_iter_chkpt = round(Int, chkpt_freq / dt)`

So the checkpoint interval is set in model time via `chkpt_freq`, but the trigger uses the iteration count; the actual spacing may differ slightly from `chkpt_freq` if `chkpt_freq / dt` is not an integer.

New runs write a `CheckpointState` (clock, grid geometry, bed, ice arrays as a dict, and basal drag `β`). Clock, grid, and ice arrays are stored as dictionaries, so extra keys are ignored and a field added later does not break reading older snapshots. A stored ice array whose shape does not match the live model is an error. Older serial files that store a `model` or `simulation` object can still be read under `BasicSpec` and `ThreadedSpec`; ice arrays are copied from the file and physics objects in the file are ignored. `MPISpec` pickup requires this snapshot format - Per-rank `Chkpt_*_Rank*` files from an earlier dev branch are not loaded.

### Checkpoint location
The directory used for both writing and reading checkpoints is resolved as follows:

1. If `chkpt_path` is set to something other than the default `"./"`, that directory is used.
2. Otherwise, if `OutputParams.output_path` is not `"./"`, checkpoints are written alongside ordinary outputs in `output_path`.
3. Otherwise, checkpoints are written under `"./"`.

For restarts, use the same `chkpt_path` and `output_path` pairing as in the run that created the checkpoint (or set `chkpt_path` explicitly to the folder that contains the `Chkpt_*.jld2` files).

Example: colocate checkpoints with field output in `outputs/`:

```julia
folder = "outputs/"
timestepping_params = TimesteppingParams(..., chkpt_freq = 1.0, chkpt_path = folder)
output_params = OutputParams(..., output_path = folder)
```

### Pickups
To continue from a checkpoint, set `niter0` to the iteration number in the filename you want to load. For example, to pick up from `Chkpt_0000001000.jld2`, use `niter0 = 1000`.

When `niter0 > 0`, `Simulation` applies `timestepping_params.dt` to the live model first, then loads `Chkpt_<niter0>.jld2` from the checkpoint directory (see above) and copies the snapshot into the `model` you passed in. The same object is kept for `BasicSpec`, `ThreadedSpec`, and `MPISpec`. Rebuild the same grid, bed, and physics as the run that wrote the file. After restore, WAVI reloads climate forcing from the last completed climate-update time on the restored clock, rebuilds surface geometry from the restored thickness, refreshes basal hydrology diagnostics without advancing water thickness, rebuilds Glen B from the restored temperature without advancing ice temperature, and rebuilds wavelets and rheological operators using that `dt`. The run then continues from iteration `niter0 + 1` up to the new `n_iter_total`.

`MPISpec` writes one checkpoint file covering the entire domain by the main rank gathering all the split data into itself before saving. For picking up from the checkpoint, it scatters into each rank's local fields, so `px`/`py` does not need to match the run that wrote the file. The same `MPISpec` snapshot can be picked up and continued with`BasicSpec`, and vice versa.

!!! note
    After a pickup, **`clock` comes from the checkpoint file** and ice arrays are copied into the **driver-built `model`**. **`TimesteppingParams` and `OutputParams` come from the arguments you pass to the new `Simulation`**. For example, a new `end_time` or `n_iter_total` controls how much further the run continues, and `output_params` controls where and how subsequent field output is written. Ensure the checkpoint file exists at the resolved checkpoint path before starting a pickup run.

    This differs from older serial checkpoints that stored a full `Model` and replaced the driver-built object on pickup. Melt, and other physics fields are no longer restored from the file; they always come from the `Model` you build in the driver.
