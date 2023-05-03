# Installation instructions

You can install the latest version of WAVI using the built-in package manager (accessed by pressing `]` in the
Julia command prompt) to add the package and instantiate/build all dependencies

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/RJArthern/WAVI.jl")

julia> Pkg.instantiate()
```

The above will install the WAVI code in the 'main' branch. To install code contained on a different branch, use the `rev` flag:


```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/RJArthern/WAVI.jl", rev="BranchName")

julia> Pkg.instantiate()
```
where `BranchName` should be replaced by the name of the branch containing the code you wish to install. 

Note that WAVI is only tested on Julia versions 1.5 and newer; stability cannot be guaranteed on newer versions!

At this time, updating should be done with care, as WAVI is under rapid development. While we take care to avoid breaking changes, they may happen during this time. If anything does break, please open an issue and let us know!