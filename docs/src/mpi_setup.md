# MPI setup and running WAVI

WAVI's distributed mode uses [`MPISpec`](./model_specifications.md#mpispec) via [MPI.jl](https://juliaparallel.org/MPI.jl/). On a new machine (laptop, workstation, or HPC cluster) you may need to configure MPI.jl before running WAVI cases via shell (using `mpiexecjl`) or cluster jobs; the built-in default however, is enough for many local and test workflows.

This page summarises the WAVI-specific workflow. For a full background on backends, and ABI compatibility, please see the [MPI.jl configuration guide](https://juliaparallel.org/MPI.jl/stable/configuration/).

!!! tip "Quick Start"
    If you just want to run WAVI across multiple cores without concerns of performance, you can skip MPI backend configuration (**skip Sections 1–2**). MPI.jl uses a built-in default library that works out of the box. For running WAVI cases via shell, jump to **[Section 3: Install `mpiexecjl`](@ref install-mpiexecjl)**; to just test that MPI is setup correctly, **[Section 4](@ref verify-installation)** is enough.

## Overview

| Step | Action                                                                                 |
| :--- | :------------------------------------------------------------------------------------- |
| 1    | Choose an MPI backend (JLL default, system MPI, or MPItrampoline).                     |
| 2    | Record your choice in this project's `LocalPreferences.toml` via `MPIPreferences`.     |
| 3    | Restart Julia and run `Pkg.instantiate()` so dependent JLLs match the MPI ABI.         |
| 4    | Verify the setup with `MPI.Get_library_version()` and a small WAVI MPI test.           |
| 5    | Install `mpiexecjl` (`MPI.install_mpiexecjl()`) for launching scripts from the shell.  |

**Tests vs Drivers:** WAVI's internal tests (such as `test/runtests_mpi_unit.jl`) spawn MPI jobs directly from Julia using `MPI.mpiexec()`. They do not require the `mpiexecjl` command. In contrast, running example drivers or your own scripts from the bash shell requires `mpiexecjl`.

**Project Configuration:** MPI configuration is set per Julia project (the directory containing `Project.toml`). Preferences in `LocalPreferences.toml` at the WAVI root apply when you use `--project=.` and to example drivers launched that way. The `test/` project has its own `Project.toml`; it usually picks up the default built-in backend unless you configure `test/` separately or root preferences propagate via Julia's load path (a root `LocalPreferences.toml` can affect `test/` too). Either default is fine for running the test suite.

## 1. Choose an MPI backend

There are three recommended options for your backend:

### Default JLL binary (easiest for a single machine)
Out of the box, MPI.jl uses `MPICH_jll`. No extra system installation is required. This is sufficient for development and to get started.

### System MPI (clusters with a site-wide install)
On HPC systems, you will often load a module (for example, `module load openmpi` or `module load intel-mpi`) and point MPI.jl at that library. You should perform this configuration in the identical environment you intend to use for jobs, ensuring the same modules and compiler wrappers are active.

```julia
using MPIPreferences
MPIPreferences.use_system_binary()   # Auto-detects the system library
# Alternatively, specify it explicitly:
# MPIPreferences.use_system_binary(abi = "OpenMPI", mpiexec = "mpiexec")
```

For more details on `LocalPreferences.toml` layout and cluster-wide preference files, see [Using a system-provided MPI backend](https://juliaparallel.org/MPI.jl/stable/configuration/#Using-a-system-provided-MPI-backend).

### MPItrampoline (recommended for safe external MPI)
[MPItrampoline](https://juliaparallel.org/MPI.jl/stable/configuration/#Using-MPItrampoline) forwards MPI calls to your site's library while keeping Julia's JLL packages ABI-compatible.

Using a system MPI directly can sometimes cause Application Binary Interface (ABI) conflicts with other pre-compiled Julia packages (such as HDF5) that expect the default MPI ABI. MPItrampoline solves this by acting as a translation layer: it guarantees compatibility with the Julia ecosystem while still routing the actual communication to your cluster's high-performance system MPI.

For installation and configuration instructions, please refer to the [MPI.jl documentation on using MPItrampoline](https://juliaparallel.org/MPI.jl/stable/configuration/#Using-MPItrampoline).

### Switching between JLL implementations
You can switch between different JLL implementations using the following command:

```julia
using MPIPreferences
MPIPreferences.use_jll_binary("OpenMPI_jll")   # or "MPICH_jll", "MPItrampoline_jll"
```

**Crucial Step:** After any change using `use_jll_binary` or `use_system_binary`, you must **quit Julia completely**. Then, reopen Julia and run `Pkg.instantiate()` in the project to ensure the changes take effect and to avoid ABI mismatch errors.

## [2. Configure the WAVI project](@id configure-wavi-project)

From the repository root (or your own project that depends on WAVI), activate the environment and set your preferences:

```julia
using Pkg
Pkg.activate(".")          # Activates the WAVI root, or your driver project
Pkg.instantiate()

using MPIPreferences
# Pick one of the options from section 1, for example:
# MPIPreferences.use_jll_binary("MPICH_jll")
# MPIPreferences.use_system_binary()
```

This process writes or updates the `LocalPreferences.toml` file in that directory. WAVI already lists `MPI` and `MPIPreferences` in its `Project.toml`, so you do not need to add them again unless you are maintaining a separate driver-only project.

## [3. Install `mpiexecjl` (shell launches only)](@id install-mpiexecjl)

You do **not** need `mpiexecjl` to run the test suite; those scripts use `MPI.mpiexec()` internally. However, you must install `mpiexecjl` if you plan to launch WAVI from a shell or a scheduler (for instance, `mpiexecjl -n 4 julia MISMIP_PLUS.jl`).

`mpiexecjl` is a Julia-aware wrapper around standard `mpiexec` or `mpirun` that ensures the correct Julia environment is used by all distributed processes. Install it once per Julia version:

```julia
using MPI
MPI.install_mpiexecjl()
```

The executable is placed in `~/.julia/bin` (assuming a standard Julia installation). You should add it to your `PATH` in Bash:

```bash
# Apply for the current shell session
export PATH="$HOME/.julia/bin:$PATH"

# Persist for future login shells (appends once to ~/.bashrc)
grep -q '\.julia/bin' ~/.bashrc 2>/dev/null || \
  echo 'export PATH="$HOME/.julia/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify the installation
which mpiexecjl
```

Use the identical `export` command in batch scripts (like Slurm) if the scheduler does not automatically load your `~/.bashrc`.

For further usage details, see [Julia wrapper for mpiexec](https://juliaparallel.org/MPI.jl/stable/usage/#Julia-wrapper-for-mpiexec).

## [4. Verify the installation](@id verify-installation)

In a **fresh** Julia session (after configuring preferences and running `instantiate`), verify that the correct library is loaded:

```julia
using MPI
MPI.Init()
println(MPI.Get_library_version())
```

From the repository root, you can test the fast MPI unit driver, which spawns two ranks internally:

```bash
julia --project=test test/runtests_mpi_unit.jl
```

Optionally, you can run one MPI test script directly using your new wrapper:

```bash
mpiexecjl -n 2 --project=test julia test/mpi_tests/test_mpispec_halo_collect.jl
```

If these commands fail with a missing `libmpi`, a wrong ABI, or launcher errors, revisit Section 1 and confirm that you fully restarted Julia after changing your preferences.

For the full list of MPI test drivers and CI behaviour, see [Unit testing: MPI tests](./test.md#mpi-tests-mpispec).

## 5. Run WAVI with `MPISpec`

Example drivers (such as [MISMIP+](./examples/mismip_plus.md)) detect MPI at startup. If they detect that `MPI.Comm_size(MPI.COMM_WORLD) > 1`, they will build an `MPISpec` and run in a distributed manner.

To run the MISMIP+ driver with four ranks:

```bash
cd example_drivers/MISMIP_PLUS
mpiexecjl -n 4 --project=../.. julia MISMIP_PLUS.jl
```

Alternatively, you can use the bundled sweep script, which tests serial, threaded, and several MPI process counts:

```bash
./run_mismip_plus.sh
```

**Grid:** `MPISpec(px, py, halo, grid)` requires the global grid at construction time. The total process count must precisely match `px * py`. For more information, see [Model specifications](./model_specifications.md). One GPU per rank: `MPISpec(..., child_architecture = GPU())` after `using CUDA`. See [GPU setup](./gpu_setup.md).

**Logging:** WAVI uses `@info` and `@debug` for progress and residuals. You can enable debug logging across distributed processes by setting the environment variable `JULIA_DEBUG=WAVI`.

## 6. HPC batch jobs

On Slurm (and similar schedulers), you must load the same MPI modules as you did in [Section 2](@ref configure-wavi-project). You then launch your script with `mpiexecjl`, deriving the task count directly from the scheduler. For example:

```bash
module load mpi/mpich-x86_64    # This is site-specific
mpiexecjl --project=../.. -n "$SLURM_NTASKS" julia MISMIP_PLUS.jl
```

A more extended, BAS-style example is available in [Running on HPC](./running_on_hpc.md#mpi-execution).

## Troubleshooting

| Symptom                                     | Likely cause                                               | Resolution                                                                                                          |
| :------------------------------------------ | :--------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------ |
| `libmpi` not found                          | System MPI is not configured, or the module is not loaded. | Run `MPIPreferences.use_system_binary()` after loading modules; check the `libmpi` path in `LocalPreferences.toml`. |
| HDF5 or other JLL errors after MPI change   | ABI mismatch across JLLs.                                  | Restart Julia, then run `Pkg.instantiate()` in every project that uses MPI.                                         |
| `mpiexecjl: command not found`              | Wrapper is not installed or not on `PATH`.                 | Run `MPI.install_mpiexecjl()`; add `~/.julia/bin` to your `PATH`.                                                   |
| Hangs or wrong rank count                   | Launcher or process layout mismatch.                       | Ensure `-n` matches `px * py`; on Slurm, use `$SLURM_NTASKS`.                                                       |

## References

- [MPI.jl — Configuration](https://juliaparallel.org/MPI.jl/stable/configuration/)
- [MPI.jl — Usage (mpiexec, writing tests)](https://juliaparallel.org/MPI.jl/stable/usage/)
- [WAVI — Model specifications (MPISpec)](./model_specifications.md)
- [WAVI — Running on HPC](./running_on_hpc.md)
- [WAVI — Unit testing (MPI tests)](./test.md#mpi-tests-mpispec)
