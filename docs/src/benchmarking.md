# Benchmarking WAVI

Compare wall time, memory, and CPU usage for [`BasicSpec`](./model_specifications.md),
[`ThreadedSpec`](./model_specifications.md), and [`MPISpec`](./model_specifications.md)
on the same example driver and grid. The harness lives under `benchmarks/` and does
not modify `src/`.

## Environment setup

Always use the **benchmarks** Julia project (not the repository root):

```bash
cd benchmarks
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

From the repository root, launch with `--project=benchmarks`:

```bash
julia --project=benchmarks benchmarks/run.jl --help
```

For MPI runs, install and configure `mpiexecjl` first, see [MPI setup](./mpi_setup.md).
On HPC, also see [Running on HPC](./running_on_hpc.md).

## CLI overview

```bash
julia --project=benchmarks benchmarks/run.jl <subcommand> [args...] [options...]
```

| Subcommand | Purpose |
|------------|---------|
| `drivers`  | List registered driver adaptors |
| `run`      | Timed benchmark with resource telemetry |
| `profile`  | On-demand CPU profile (`@profile`) |
| `plot`     | Overlay RSS/CPU time series from one or more runs |
| `mpi-overhead` | Calculate pure MPI setup overhead from JSON results |

Each subcommand supports `--help`.

## Running timed benchmarks

```bash
# BasicSpec (serial baseline)
julia --project=benchmarks -t 1 benchmarks/run.jl run basic mismip_plus

# ThreadedSpec (shared-memory Schwarz; -t ≥ ngridsx × ngridsy)
julia --project=benchmarks -t 4 benchmarks/run.jl run threaded mismip_plus --ngridsx 2 --ngridsy 2

# MPISpec (one Julia thread per rank; prefer a 1D layout on narrow domains like MISMIP+)
mpiexecjl --project=benchmarks -n 4 julia -t 1 benchmarks/run.jl run mpi mismip_plus --px 4 --py 1
```

MISMIP+ has a thin `y` extent, so `--py > 1` leaves too few core cells after halo for efficient PoU.
Default CLI layout is `N×1` (`--px 0` uses the MPI world size with `--py 1`).

Useful `run` options: `--niterations`, `--overlap`, `--sample-interval`,
`--no-plots`, `--warmup`. See `julia --project=benchmarks benchmarks/run.jl run --help`.

Results are written under `benchmarks/output/<driver_name>/benchmark_<id>_<timestamp>/`
(`benchmark_results.json`, `resource_timeseries.csv`, a copy of the `<driver_name>.jl` adaptor, and optional plots).

## Native code caching (This is planned behaviour - implementation was removed in latest changes due to CI failure).

WAVI uses [PrecompileTools.jl](https://github.com/JuliaLang/PrecompileTools.jl) so that installing or
`Pkg.precompile()`ing the package caches native code for a small toy `BasicSpec` /
`ThreadedSpec` run (and, when MPI allows it during precompile, a single-rank `MPISpec`).

This reduces cold start for core solvers. It does **not** remove all first-use compile on multi-rank
MPI jobs: methods specialised only for multi-process layouts, and paths the toy workload never hits,
can still JIT once. After a successful precompile, a harness `--warmup` (or any short first run in
the same process) is the practical way to finish warming MPI ranks before timed iterations.

To rebuild the cache after pulling large `src/` changes:

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Set `WAVI_PRECOMPILE_MPI=0` if package precompile hangs on MPI init in your environment.

## Drivers

Drivers are thin adaptors in `benchmarks/drivers/`. Current example adaptors include:

| Driver        | Example wrapped                         |
|---------------|-----------------------------------------|
| `mismip_plus` | `example_drivers/MISMIP_PLUS/`          |
| `iceberg`     | `example_drivers/Iceberg/Iceberg.jl`    |

List registered names with:

```bash
julia --project=benchmarks benchmarks/run.jl drivers
```

To add a driver: create `benchmarks/drivers/<name>.jl` with module name `<name>`,
point at an example driver, and export `driver_run`, `driver_grid`, and
`driver_plot_vars`. No modification of the benchmarking code is required, it is
picked up automatically.

Since these cases have very small domain sizes, they are only representative of the approach, and you should utilise larger cases to do full benchmarks.

## BLAS and threads

The harness pins BLAS to one thread so OpenBLAS/MKL do not oversubscribe cores
alongside Julia threads or MPI ranks. Prefer launching with the `-t` counts above,
and set `OPENBLAS_NUM_THREADS=1` / `MKL_NUM_THREADS=1` in the environment when
running on HPC.

| Spec         | Julia `-t`             | Notes                               |
|--------------|------------------------|-------------------------------------|
| BasicSpec    | `1`                    | Serial baseline                     |
| ThreadedSpec | `≥ ngridsx × ngridsy`  | Single node only                    |
| MPISpec      | `1` per rank           | `px × py` must equal MPI world size |

## Telemetry

During `run`, an external sampler records RSS and CPU over time (so sampling still
works under `julia -t 1`). Each run directory includes:

- `benchmark_results.json` — wall time, peak RSS, allocations, `reference_cores`
  - Includes `setup_time_seconds` and `solve_time_seconds` to isolate solver math performance.
  - Includes `compilation_and_overhead_time_seconds` which captures all remaining time (JIT compilation, MPI launch, file I/O).
- `resource_timeseries.csv` — elapsed time vs RSS and CPU (including per-processor breakdowns for MPI runs)

Compare series visually:

```bash
julia --project=benchmarks benchmarks/run.jl plot \
  benchmarks/output/mismip_plus/benchmark_basic_* \
  benchmarks/output/mismip_plus/benchmark_threaded_* \
  benchmarks/output/mismip_plus/benchmark_mpi_*
```

For MPI runs, the plot will break down memory usage to show Rank 0 vs the average of other processors to help identify bottlenecks.

### Calculate MPI Setup Overhead

To isolate the exact cost of establishing the MPI topology and halo exchange buffers compared to a standard BasicSpec run, you can compare the `setup_time_seconds` of a pair of serial and MPI runs (with identical configuration):

```bash
julia --project=benchmarks benchmarks/run.jl mpi-overhead \
  benchmarks/output/mismip_plus/benchmark_basic_20260814_120000 \
  benchmarks/output/mismip_plus/benchmark_mpi_20260814_120500
```
