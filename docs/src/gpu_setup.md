# GPU setup and running WAVI

WAVI's single-GPU mode uses [`GPUSpec`](./model_specifications.md#gpuspec) through [CUDA.jl](https://cuda.juliagpu.org/). CUDA is a *weak dependency* which means that installing WAVI does not install CUDA, and the core package still runs on CPU if CUDA is not present. To use NVIDIA GPUs, you need to add CUDA to your project, load it, then pass `GPUSpec()` as `model.spec`.

Currently, only NVIDIA devices are supported. The kernels are written using KernelAbstractions.jl with the aim of adding other vendors (AMD, Intel) later without having to rewrite the core code. This relies on the extension approach added since Julia 1.9.

!!! tip "Quick start"
    On a machine with a working NVIDIA driver and GPU:

    ```julia
    using Pkg
    Pkg.add("CUDA")
    using CUDA
    CUDA.functional()   # should return true
    using WAVI

    grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
    bed = zeros(grid.nx, grid.ny)
    # Either
    model = Model(grid, bed, GPUSpec())
    # or
    model = Model(;grid=grid, bed_elevation=bed, spec=GPUSpec())
    ```

    If `GPUSpec()` results in an error, the vendor package is not loaded or no functional GPU was found. See [Section 3](@ref verify-cuda).

## Overview

| Step | Action |
| :--- | :----- |
| 1    | Add `CUDA` to the Julia project that runs your driver (not to WAVI itself). |
| 2    | Restart Julia, `using CUDA`, and check `CUDA.functional()`. |
| 3    | Build the model with `spec = GPUSpec()`. |
| 4    | On a cluster, request a GPU from the scheduler. |

**WAVI does not list CUDA in `[deps]`.** That keeps CPU CI and machines without NVIDIA hardware working. The CUDA extension (`WAVICUDAExt`) loads only after `using CUDA`.

## [1. Add CUDA.jl to your project](@id add-cuda-jl)

From the project that contains your driver (`Project.toml` in that directory):

```julia
using Pkg
Pkg.activate(".")
Pkg.add("CUDA")
Pkg.instantiate()
```

Do **not** add CUDA to WAVI's own `Project.toml`.

After adding the package, **quit Julia and start a new session** so the WAVI CUDA extension can load.

## [2. Load CUDA before WAVI GPU types](@id load-cuda)

`GPU()` and a functional `gpu_device()` exist only after CUDA is loaded *and* a device is present:

```julia
using CUDA
using WAVI

CUDA.functional()          # true if a device can run kernels
CUDA.versioninfo()         # driver, runtime, and device summary
gpu_device()               # WAVI GPU architecture, or CPU() with a warning
```

Order matters: `using CUDA` must happen in the session before `GPUSpec()`. A convenience path is to put `using CUDA` in your driver or in `~/.julia/config/startup.jl` (or, wherever this config is set-up) on GPU systems.

## [3. Verify the installation](@id verify-cuda)

In a **fresh** Julia session, with your driver project activated, first check that CUDA can run a kernel:

```julia
using CUDA
x = CUDA.ones(4)
x .+= 1
Array(x)   # [2, 2, 2, 2]
```

Then check that WAVI can see the GPU and put fields on the device:

```julia
using CUDA
using WAVI

@assert CUDA.functional()
spec = GPUSpec()
architecture(spec)         # GPU{...}, not CPU()

grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
bed = zeros(grid.nx, grid.ny)
model = Model(grid, bed, spec)
architecture(model)        # GPU
typeof(model.fields.gh.h)  # CuArray, not Matrix
```

If `GPUSpec()` throws `ArgumentError` about needing a GPU, either CUDA is not loaded, `CUDA.functional()` is false, or you are on a node with no device.

## [4. Run WAVI with GPUSpec](@id run-gpuspec)

```julia
using CUDA
using WAVI

grid = Grid(nx = 256, ny = 128, dx = 1.0e3, dy = 1.0e3)
bed = zeros(grid.nx, grid.ny)   # still a host Array at construction
model = Model(grid, bed, GPUSpec())

update_state!(model)
```

`GPUSpec()` calls `gpu_device(; force = true)`, so it errors instead of silently falling back to CPU. This is intentional so that a SLURM job which is expected to be on GPU doesn't unexpectedly and silently run on a CPU.

## [5. HPC batch jobs](@id hpc-batch)

`GPUSpec` is one process on one GPU. Several GPUs on one node (or one GPU per node) use `MPISpec` with `child_architecture = GPU()`. Halo exchange and Partition of Unity strips copy through host MPI buffers (this should be optimised in the future, currently not using CUDA-aware MPI). CUDA.jl includes the Nvidia toolkit, so, do not run a `module load cuda` on a HPC.

### One GPU (`GPUSpec`)

Pin the job to a GPU partition and start Julia as usual (no `mpiexecjl`). The following is the set-up for a GPU node on BAS HPC.

Example Slurm fragment (names of partitions and modules are site-specific):

```bash
#SBATCH --partition=gpu
#SBATCH -N 1
#SBATCH --nodelist=bsl-node-s22
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00

julia --project=. driver.jl
```

In `driver.jl`, load CUDA before constructing the model, as in [Section 4](@ref run-gpuspec).

`GPUSpec` uses one GPU. If the job can see several devices (for example `--gres=gpu:2` or a whole node), set `CUDA_VISIBLE_DEVICES` to a single index. Do not override it on a `--gres=gpu:1` job: Slurm has already pointed the process at the allocated GPU.

### Multiple GPUs (`MPISpec` + `GPU()`)

This approach can be used where you want to:
* Run on one system with multiple GPUs.
* Run on multiple systems with one or more GPUs per system.

For this, ask the scheduler for one GPU per MPI rank, and leave `CUDA_VISIBLE_DEVICES` unset. WAVI then gives each rank a GPU on that machine: the first rank on the node uses GPU 0, the second uses GPU 1, and so on. If you launch more ranks than GPUs, they share cards and WAVI prints a warning. That warning is not shown when Slurm has already given each rank its own GPU (for example `--gpus-per-task=1`).

```bash
#SBATCH --partition=gpu
#SBATCH -N 1
#SBATCH --nodelist=bsl-node-s22
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00

mpiexecjl -n 4 --project=. julia driver.jl
```

```julia
using CUDA
using WAVI
spec = MPISpec(4, 1, 2, grid; child_architecture = GPU())
model = Model(grid, bed, spec)
```

Do not combine `child_architecture = GPU()` with `local_spec = ThreadedSpec(...)`.

## Troubleshooting

| Symptom | Likely cause | Resolution |
| :------ | :----------- | :--------- |
| `GPUSpec()` throws “needs a GPU architecture” | CUDA not loaded, or no functional device. | `using CUDA` in this session; check `CUDA.functional()` and that you have a CUDA device. |
| `CUDA.functional()` is `false` | Driver, module, or visible devices. | GPU node, site CUDA module, `CUDA_VISIBLE_DEVICES`. |
| `GPU()` is a `MethodError` | CUDA extension did not load. | Add `CUDA` to *this* project; restart Julia; `using CUDA` then `using WAVI`. |
| Fields are still `Array` after `Model(...)` | Spec is not `GPUSpec`, or `MPISpec` was built without `child_architecture = GPU()`. | Pass `GPUSpec()` as the third argument to `Model`, or `MPISpec(..., child_architecture = GPU())`. Default `Model(grid, bed)` is CPU `BasicSpec`. |
| Job hangs or uses the wrong GPU on a multi-GPU node | Several devices visible to one `GPUSpec` process. | Set `CUDA_VISIBLE_DEVICES` to one index for `GPUSpec` only. For MPI+GPU, leave it unset so ranks pin by node-local rank. |
| `ThreadedSpec` with `child_architecture = GPU()` | Invalid combination. | Use `local_spec = nothing` (the default) on GPU ranks. |

## References

- [CUDA.jl](https://cuda.juliagpu.org/)
- [WAVI — Model specifications (GPUSpec)](./model_specifications.md)
- [WAVI — Running on HPC](./running_on_hpc.md)
- [WAVI — MPI setup](./mpi_setup.md)
