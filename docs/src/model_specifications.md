# WAVI Model Specifications

## Overview

A late addition to WAVI is the idea of specifications, which guide the setup of the library. 

Semantically, we can think of the components of WAVI this way: 

* `Model` brings together all the elements of compute and data into a singular operable entity
* `Simulation` determines how the model will be used to produce results over time
* `Specs` determines how the model and simulation work together to solve the ice sheet system

Specifications heavily leverage multiple dispatch to re-implement portions of WAVI to introduce new structural architectures, computational approaches and data operations.

## [Choosing how to run WAVI](@id choosing-how-to-run-wavi)

`spec` is how you tell WAVI *where* the local solve runs and *whether* the domain is split for parallel processing. The type names in code are `BasicSpec`, `ThreadedSpec`, `MPISpec`, and `GPUSpec`.

| Option | Spec | How you launch | Memory | Use when |
| :--- | :--- | :--- | :--- | :--- |
| Serial CPU | `BasicSpec()` (default) | `julia --project=. driver.jl` | one process | development, small domains |
| CPU threads | `BasicSpec()` | `julia -t N --project=. driver.jl` | one process | one node; KernelAbstractions.jl uses Julia threads |
| Shared-memory Schwarz | `ThreadedSpec(...)` | `julia -t N` with `N` at least `ngridsx * ngridsy` | one process | overlapping subdomain preconditioner; still one node's RAM |
| Distributed CPU | `MPISpec(px, py, halo, grid)` | `mpiexecjl -n px*py ...` | per rank | several nodes, or several processes on one node |
| Single GPU | `GPUSpec()` | `julia --project=. driver.jl` after `using CUDA` | one GPU | NVIDIA GPU; see [GPU setup](./gpu_setup.md) |
| One GPU per MPI rank | `MPISpec(..., child_architecture = GPU())` | `mpiexecjl -n N ...` after `using CUDA` | per rank, one GPU | several GPUs on one node, or one GPU per node |

`BasicSpec` plus `julia -t N` is **not** the same as `ThreadedSpec`. Threads on `BasicSpec` parallelise loops on the whole grid. `ThreadedSpec` builds overlapping Schwarz subdomains and solves each on a Julia thread which can result in slightly varying results.

`GPUSpec` is one GPU in one process. Several GPUs use `MPISpec` with `child_architecture = GPU()` (one rank, one GPU). Halo exchange copies through host MPI buffers; this is not CUDA-aware MPI. Do not pass `local_spec = ThreadedSpec(...)` with a GPU child. Use [GPU setup](./gpu_setup.md) for CUDA.jl, and [MPI setup](./mpi_setup.md) for `MPISpec`.

### Implementation notes
* **Grid balancing:** if the global grid is not divisible by the number of subdomains or MPI ranks, remainder cells are spread over the first few subdomains (`ThreadedSpec` and `MPISpec`). You do not have to pick a thread or rank count that divides the grid exactly.
* Examples in these docs default to `BasicSpec`, so a model is not distributed unless you pass another spec.

## Setting up specifications

Pass a [`Grid`](@ref) and bed into [`Model`](@ref). With no `spec`, a `BasicSpec` is used.

```julia
model = Model(grid, bed)   # BasicSpec, CPU

spec = ThreadedSpec(
    ngridsx = 16,
    ngridsy = 2,
    overlap = 1,
    niterations = 1,
)
model = Model(grid, bed, spec)

# px, py, halo, grid  (px * py must equal the MPI world size)
spec = MPISpec(16, 2, 1, grid)
model = Model(grid, bed, spec)

using CUDA
model = Model(grid, bed, GPUSpec())

spec = MPISpec(4, 1, 2, grid; child_architecture = GPU())
model = Model(grid, bed, spec)
```

Pass `spec` as the third argument, or use `Model(; grid, bed_elevation, spec=...)`. `Model(grid, bed; spec=...)` will not work: the two-argument method always inserts `BasicSpec()`.

The grid is required when you construct `MPISpec`, because the split happens then. `ThreadedSpec` and `GPUSpec` do not need the grid at spec construction.

### Using multiple CPU threads (`BasicSpec`)

```bash
julia -t 8 --project=. driver.jl
```

No change in specification, just the default `BasicSpec`. KernelAbstractions uses those Julia threads to parallelise across loops.

### Threaded Schwarz (`ThreadedSpec`)

Use `ThreadedSpec` as above, and launch Julia with enough threads for the subdomain grid (`ngridsx * ngridsy`).

### MPI execution

Configure MPI once per machine (MPIPreferences, `mpiexecjl`, checks). See [MPI setup](./mpi_setup.md) and the [MPI.jl configuration guide](https://juliaparallel.org/MPI.jl/stable/configuration/).

```bash
mpiexecjl -n <num> --project=<path-to-project> julia <path-to-driver.jl> <driver args...>
```

Example (MISMIP+ driver, four ranks):

```bash
mpiexecjl -n 4 --project=../.. julia example_drivers/MISMIP_PLUS/MISMIP_PLUS.jl
```

### GPU execution

Add CUDA.jl to the driver project and load it. One GPU: `GPUSpec()`. Several GPUs: one MPI rank per GPU with `MPISpec(..., child_architecture = GPU())`. Step-by-step: [GPU setup](./gpu_setup.md).

```bash
julia --project=. driver.jl    # In the driver: using CUDA; GPUSpec()
mpiexecjl -n 4 --project=. julia driver.jl    # 4 ranks in this case means 4 GPUs; In the driver: MPISpec(..., child_architecture = GPU())
```

## Domain Decomposition

### ThreadedSpec

The threaded specification utilises an iterative Schwarz approach to domain decomposition, producing individual domains prior to execution of the velocity solve processing. 

For this specification only two methods are overridden:

`update_preconditioner!` creates a decomposed grid across individual threads, allowing these threads to `precondition!` individually within that subdomain.

`precondition!` then handles exchange of velocities takes place during each iteration of the velocity solve. The state is updated following a transfer of velocities from global to the "local" grid, upon which the local `update_state!` is used to update the local model prior to transferring back to the global grid using a [partition of unity](https://en.wikipedia.org/wiki/Partition_of_unity).

The memory space is still limited to one process, with the potential to leverage threading optimisations within a single physical processor. _Therefore this might increase parallelisation of the model computations, __but it will remain memory bound__._

### GPUSpec

`GPUSpec` does not split the domain. One process owns the whole grid; dense arrays live on one NVIDIA GPU. Some aspects still live and are transferred across from host (CPU) to device (GPU) which should be optimised for in the future. Colour-list Gauss-Seidel and Kronecker geometry (`cent`, `∂x`, `∂y`) still run on the host and copy the fields they need.

Use [`GPUSpec()`](./gpu_setup.md) after `using CUDA`. Several GPUs use [`MPISpec`](#mpispec) with `child_architecture = GPU()`.

### MPISpec

The MPI specification is used to decompose the global domain of the `Model`. Unlike the threaded specification however, this decomposition takes place at time of creation. The grid provided is split and the root rank takes the upper-left node in the topology, [which is clearly explained here](https://hpc-tutorials.llnl.gov/mpi/virtual_topologies/).

The list of methods overridden is not explained in detail here, there are numerous areas of the model that require alteration:

* The underlying `Model` constructor is implemented such that each node contains it's own localised portion of the global domain.
* Pass `child_architecture = GPU()` after `using CUDA` to put each rank's local fields on one GPU. Rank 0's `global_fields` (output gather) stay on the host. Halo exchange is host-staged (not CUDA-aware MPI). `local_spec = ThreadedSpec(...)` is rejected on GPU ranks.

!!! tip "Process Layout Geometry"
    The geometric layout of the MPI process grid (`px` × `py`) heavily impacts performance and solver convergence. For domains with high aspect ratios (like MISMIP+ which is long and narrow), **favour 1D process layouts** (e.g., `4x1` instead of `2x2`).

    A `2x2` layout on a narrow grid creates subdomains with very few "core" cells compared to the halo size, which drastically increases the relative overhead of Schwarz Partition-of-Unity (`pou`) communication and slows down convergence.

* Velocity solves use overlapping **Schwarz iterations** (multiple local subdomain solves per step). The subdomain boundary treatment is controlled by the `pou` parameter:
    * **`pou=true` (default)**: Uses a **Partition of Unity** approach. Subdomain velocities are smoothly blended together in overlap zones, ensuring physical continuity. This requires a `halo` size of at least `2`.
    * **`pou=false`**: Uses **Restricted Additive Schwarz (RAS)**. Velocities are exchanged directly between neighboring subdomains without blending. This is computationally simpler and allows a smaller `halo` size of `1`.
* **Thickness Halos**: Ice thickness `h` is synchronised across subdomains at the end of each timestep.
* **Global Residual Check**: Picard convergence is verified globally across all ranks. Each process computes the squared norm of its local residual on "core" unknowns only (excluding halo overlaps) to prevent double-counting. These are then aggregated via `MPI.Allreduce` to determine the global relative residual.
* `Simulation` operations are updated to collect data as required for outputs using a deferred evaluation mechanism, required to ensure outputs are taken from the global grid and rendered from the root process only.

The topology of an individual ranks grid (in this case rank zero) is also relevant, with velocities at the neighboured-edges fixed and then transferred. The following diagram illustrates the situation:

```@raw html
<center><img src="../assets/mpispec_exchange.png" alt="" title="" height="700" /></center>
```

#### Distributed TODO (MPISpec)

There are elements of MPISpec still under development at time of writing: 

* [Outputs only can be taken from `model.global_fields`](https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues/111) - therefore melt rates, parameters and other fields are not accessible as outputs (though `model.mpi_rank` is now available).
* Non-2D fields (e.g., 3D temperature) are not yet accessible as distributed outputs.

## Specifications API

Please refer to the [API documentation for the structure / constructor definitions](API/specifications.md).