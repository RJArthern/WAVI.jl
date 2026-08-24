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

Currently, only running on one GPU is supported while MPI + CUDA is planned. Pin the job to a GPU partition and start Julia as usual (no `mpiexec` for `GPUSpec`). The following is the set-up for a GPU node on BAS HPC.

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

## Troubleshooting

| Symptom | Likely cause | Resolution |
| :------ | :----------- | :--------- |
| `GPUSpec()` throws “needs a GPU architecture” | CUDA not loaded, or no functional device. | `using CUDA` in this session; check `CUDA.functional()` and that you have a CUDA device. |
| `CUDA.functional()` is `false` | Driver, module, or visible devices. | GPU node, site CUDA module, `CUDA_VISIBLE_DEVICES`. |
| `GPU()` is a `MethodError` | CUDA extension did not load. | Add `CUDA` to *this* project; restart Julia; `using CUDA` then `using WAVI`. |
| Fields are still `Array` after `Model(...)` | Spec is not `GPUSpec`. | Pass `GPUSpec()` as the third argument to `Model`. Default `Model(grid, bed)` is CPU `BasicSpec`. |
| Job hangs or uses the wrong GPU on a multi-GPU node | Several devices visible. | Set `CUDA_VISIBLE_DEVICES` to one index only if the job can see more than one GPU. |

## References

- [CUDA.jl](https://cuda.juliagpu.org/)
- [WAVI — Model specifications (GPUSpec)](./model_specifications.md)
- [WAVI — Running on HPC](./running_on_hpc.md)
- [WAVI — MPI setup](./mpi_setup.md)
