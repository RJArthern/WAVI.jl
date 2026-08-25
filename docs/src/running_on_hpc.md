# Running on HPC Clusters

WAVI.jl is designed to be highly scalable and can be run efficiently on High-Performance Computing (HPC) clusters. Below are examples of how to submit jobs using the SLURM workload manager for serial, multi-threaded, MPI, and GPU configurations.

These examples are based on the `MISMIP_PLUS` driver provided in the `example_drivers/MISMIP_PLUS` directory for running on the British Antarctic Survey's HPC.

## Serial Execution

For smaller test runs or non-parallel configurations, you can run WAVI.jl in serial. 

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.serial.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```

## Multi-threaded Execution

To take advantage of shared-memory parallelism on a single node, you can run Julia with multiple threads. Set the `--cpus-per-task` in your SLURM script and pass this value to Julia using the `-t` flag.

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.threaded.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```

## MPI Execution

For large-scale simulations spanning multiple nodes or making use of distributed memory, WAVI.jl supports MPI via [`MPISpec`](./model_specifications.md#mpispec).

!!! note "Before your first MPI job"
    Configure MPI.jl, install `mpiexecjl`, and verify with the MPI unit tests. Step-by-step instructions are in [MPI setup](./mpi_setup.md).

On the cluster, load the site MPI module if that was your approach to setting up mpi (e.g. `mpich`) in the batch script, then launch with `mpiexecjl` and the scheduler task count (`$SLURM_NTASKS` below).

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.mpi.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```

## GPU execution

See the [GPU setup](./gpu_setup.md) section to understand how to request GPUs, load vendor packages (CUDA, ... more in the future), run on individual GPUs using `GPUSpec`, and a combination of `MPISpec` and `GPUSpec` to run on multiple GPUs.
