# WAVI Benchmarks

This directory contains the benchmarking suite for WAVI.jl, allowing you to trace performance improvements and profile CPU/Memory usage.

## Running on BAS HPC (SLURM)

### CPU scaling (`run_scaling.sh`)

The `run_scaling.sh` script is for CPU-only multi-core scaling sweeps (BasicSpec thread counts and MPISpec ranks) on the British Antarctic Survey (BAS) HPC. It is configured for a compute node with 36 physical cores on the `medium` partition.

### Usage

You can run the sweep directly (if already on compute node):
```bash
bash benchmarks/run_scaling.sh "my_custom_tag"
```

Alternatively, you can submit it directly to the SLURM queue from a login node:
```bash
sbatch benchmarks/run_scaling.sh "my_custom_tag" "mismip_plus"
```

*(If you omit the arguments, they default to the values in the script)*

To run a simple single-core `BasicSpec` (serial) baseline locally or interactively:

```bash
julia -t 1 --project=benchmarks benchmarks/run.jl run basic ismip7_16km_synthetic --tag "baseline_sparse"
```

To run a simple single-core baseline on a compute node directly from a login node:

```bash
srun -A short -p short -N 1 -n 1 --time=00:30:00 julia -t 1 --project=benchmarks benchmarks/run.jl run basic ismip7_16km_synthetic --tag "baseline_sparse"
```

### GPU trial run (`run_gpu.sh`)

`run_gpu.sh` is a one-GPU `GPUSpec` job (not a thread or MPI sweep). The script defaults to `bsl-node-s22` (which currently has 8 x V100S). CUDA.jl vendors its own toolkit (So, you must not do a `module load cuda`), the NVIDIA driver on the GPU node is enough. It requests `--gres=gpu:1` and does not use `mpiexec`. CUDA.jl is added to the benchmarks project on first run if it is missing. `sbatch --nodelist=` overrides the default node.

```bash
sbatch benchmarks/run_gpu.sh "ka.jl_gpuspec" "ismip7_16km_synthetic"
sbatch --nodelist=bsl-node-s20 benchmarks/run_gpu.sh "ka.jl_gpuspec"   # Currently has 2 A40 GPUs
sbatch --nodelist=bsl-node-s21 benchmarks/run_gpu.sh "ka.jl_gpuspec"   # Currently has 2 T4 GPUs
```

If you are already on the GPU node:

```bash
bash benchmarks/run_gpu.sh "ka.jl_gpuspec"
```

Or, without the script:

```bash
julia --project=benchmarks -t 1 benchmarks/run.jl run gpu ismip7_16km_synthetic --tag "ka.jl_gpuspec_gpu1"
```

### Outputs

The CPU script runs a scaling sweep (1, 2, 4, 8, 16, and 36 cores) using the `ismip7_16km_synthetic` driver. The GPU script runs the same driver once with `GPUSpec`. Telemetry is saved under `benchmarks/output/<driver>/`.

Each run generates:
- `resource_timeseries.csv`: CPU and memory (RSS) usage sampled at 0.5s intervals.
- `benchmark_results.json`: Metadata about the run, including the exact Git commit hash and the `--tag` (e.g., `baseline_sparse_36core`).
- Standard NetCDF outputs and basic diagnostic plots.
