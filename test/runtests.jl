group = "all" 

if group == "fields" || group == "all"
    println(1)
    include("test_fields.jl")
end

if group == "grids" || group == "all"
    include("test_grids.jl")
end

if group == "melt" || group == "all"
    include("test_melt.jl")
end

if group == "models" || group == "all"
    include("test_models.jl")
end

if group == "outputs" || group == "all"
    include("test_outputting.jl")
end

if group == "simulations" || group == "all"
    include("test_simulations.jl")
end

if group == "timesteppingparams" || group == "all"
    include("test_timesteppingparams.jl")
end

if group == "utils" || group == "all"
    include("test_utils.jl")
end

if group == "kronecker" || group == "all"
    include("test_kronecker.jl")
end

if group == "version_update" || group == "all"
    include(joinpath("version_update_test_verification","test_version_updates.jl"))
end

if group == "shared_memory_parallelism" || group == "all"
    include(joinpath("domain_decomposition_tests","test_domains_shared_memory.jl"))
end

if group == "verification" || group == "all"
    include("verification_tests.jl")
end