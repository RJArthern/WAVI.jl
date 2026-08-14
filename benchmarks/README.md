# WAVI Benchmarks

This directory contains the benchmarking suite for WAVI.jl, allowing you to trace performance improvements and profile CPU/Memory usage.

## Running on BAS HPC (SLURM)

The `run_scaling.sh` script is tailored for running multi-core scaling sweeps on the British Antarctic Survey (BAS) HPC using SLURM. It is configured for the compute node with 36 physical cores and using the `medium` partition.

### Usage

You can run the sweep directly (if already on compute node):
```bash
bash benchmarks/run_scaling.sh "my_custom_tag"
```

Alternatively, you can submit it directly to the SLURM queue from a login node:
```bash
sbatch benchmarks/run_scaling.sh "my_custom_tag"
```

*(If you omit the tag argument, it defaults to `"baseline_sparse"`)*

To run a simple single-core `BasicSpec` (serial) baseline locally or interactively:

```bash
julia -t 1 --project=benchmarks benchmarks/run.jl run basic ismip7_16km_synthetic --tag "baseline_sparse"
```

To run a simple single-core baseline on a compute node directly from a login node:

```bash
srun -A short -p short -N 1 -n 1 --time=00:30:00 julia -t 1 --project=benchmarks benchmarks/run.jl run basic ismip7_16km_synthetic --tag "baseline_sparse"
```

### Outputs

The script runs a scaling sweep (1, 2, 4, 8, 16, and 36 cores) using the `ismip7_16km_synthetic` driver. All telemetry data and outputs are saved to `benchmarks/output/ismip7_16km_synthetic/`.

Each run generates:
- `resource_timeseries.csv`: CPU and memory (RSS) usage sampled at 0.5s intervals.
- `benchmark_results.json`: Metadata about the run, including the exact Git commit hash and the `--tag` (e.g., `baseline_sparse_36core`).
- Standard NetCDF outputs and basic diagnostic plots.
