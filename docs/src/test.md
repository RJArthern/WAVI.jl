# WAVI.jl tests
WAVI.jl includes a comprehensive suite of unit tests, which help to prevent bugs, ensure code quality, and prevent breaking changes. Scripts to run these tests are held in the `/test` folder. These tests utilize the [Julia unit testing architechture](https://docs.julialang.org/en/v1/stdlib/Test/), which provides a streamlined environment for testing.

WAVI.jl includes 10 individual testing scripts, `/test/test_XXXXX.jl`. These roughly correspond to each level of the [data structure](./API/overview.md) heirarchy, with a couple of extras. They are
- `/test/test_fields.jl`: test the [WAVI.jl fields structures](./API/fields.md).
- `/test/test_grids.jl`: test the [WAVI.jl grid structures](./API/fields.md).
- `/test/test_kronecker.jl`: test the Kronecker product operations used in WAVI.jl.
- `/test/test_melt.jl`: test the [WAVI.jl basal melt rate parametrizations](./API/melt_rate_models.md).
- `/test/test_models.jl`: test the [WAVI.jl model structures](./API/model.md).
- `/test/test_outputting.jl`: test the [WAVI.jl output writing structures](./API/output_writing.md).
- `/test/test_simulations.jl`: test the [WAVI.jl simulation structures](./API/simulations.md).
- `/test/test_timesteppingparams.jl`: test the [WAVI.jl timestepping parameters structures](./API/timestepping_params.md).
- `/test/test_utils.jl`: test the utility functions employed by WAVI.jl.
- `/test/verification_tests.jl`: verification tests for the WAVI.jl model. Note that these tests require setup and running of a WAVI.jl simulation and therefore take a while!   

Individual unit tests can be run by running the corresponding Julia script (note that WAVI.jl must be installed)
```julia
julia> include('./test/test_XXXXX.jl')
```
The output will appear as below (in this case for `/test/test_fields.jl`). The test summary lists the total number of tests, alongside the number which pass and fail. 

```@raw html
<center><img src="https://raw.githubusercontent.com/RJArthern/WAVI.jl/docs-reconcile/docs/src/assets/test_output.png" alt="" title="" width="800"  /></center>
```

The full suite of tests can be run using the `runtests.jl` script. 
```julia
julia> include('./test/runtests.jl')
```

## Continuous Integration