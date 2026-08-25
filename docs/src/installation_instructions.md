# Installation instructions

!!! details "Installing Julia"
    We recommend [juliaup](https://github.com/JuliaLang/juliaup) to install and switch Julia versions on your machine. It is the official cross-platform Julia version manager (replacing manual downloads or distribution packages for day-to-day use).

    * **Install juliaup** (see the [juliaup README](https://github.com/JuliaLang/juliaup#installation) for Windows, macOS, and Linux options). On Linux/macOS, a typical one-line install is:
      ```bash
      curl -fsSL https://install.julialang.org | sh
      ```
      Restart your shell so `julia` and `juliaup` are on your `PATH` (often `~/.juliaup/bin`).

    * **Add Julia versions** used by WAVI development and CI:
      ```bash
      juliaup add 1.11
      juliaup add 1.12   # optional: matches structural-test matrix in CI
      juliaup default 1.11
      ```
      Check what is installed:
      ```bash
      juliaup status
      julia --version
      ```

    * **Use a specific version** for one command without changing the default:
      ```bash
      julia +1.11 --project -e 'using Pkg; Pkg.instantiate()'
      ```
      When working from a git clone, `cd` into the repo and use `--project` (or `--project=test` for tests) so dependencies resolve against that environment.

    WAVI is tested in CI on Julia **1.11** and **1.12** (structural tests); MPI unit tests run on **1.12**.

## Overview

You can install the latest version of WAVI using the built-in package manager to add the package and instantiate/build all dependencies:

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/WAVI-ice-sheet-model/WAVI.jl")

julia> Pkg.instantiate()
```

or, alternatively, you can use the shortcut `]` to trigger the Pkg REPL and run these simplified commands:

```julia
julia> ]

(@v1.12) pkg> add https://github.com/WAVI-ice-sheet-model/WAVI.jl

(@v1.12) pkg> instantiate
```

## Branch selection

The above will install the WAVI code in the 'main' branch. To install code contained on a different branch, use the `rev` flag:


```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/WAVI-ice-sheet-model/WAVI.jl", rev="BranchName")

julia> Pkg.instantiate()
```
where `BranchName` should be replaced by the name of the branch containing the code you wish to install. 

Note that WAVI is tested in CI on Julia **1.11** and **1.12**; other versions may work but are not guaranteed on every branch.

At this time, updating should be done with care, as WAVI is under rapid development. While we take care to avoid breaking changes, they may happen during this time. If anything does break, please open an issue and let us know!

## MPI

Distributed runs use [MPI.jl](https://juliaparallel.org/MPI.jl/) and [`MPISpec`](./model_specifications.md#mpispec). On a new machine you must [configure MPI](./mpi_setup.md) before running WAVI on multiple processes.

## GPU

Single-GPU runs use [CUDA.jl](https://cuda.juliagpu.org/) and [`GPUSpec`](./model_specifications.md#gpuspec). CUDA is not a WAVI dependency: add it to your driver project, then [set up the GPU](./gpu_setup.md).
