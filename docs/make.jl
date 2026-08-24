using Documenter
using DocumenterCitations
using Literate 
using Plots
using Revise
using WAVI

# Track changes in WAVI source code
# Handle need for re-compiling instead of using previous cache
Revise.track(WAVI)

bib_filepath = joinpath(dirname(@__FILE__), "wavi.bib")
bib = CitationBibliography(bib_filepath)

@info "Running makedocs()"

ENV["GKSwstype"] = "100"

#####
#### Organize page hierarchies
#####

#use these examples if you aren't generating on the fly
example_pages = [
    "First example: one-dimensional planar flow"    => "examples/planar_1D_flow.md"
    "MISMIP+ part one: grounding lines on overdeepened bed" => "examples/mismip_plus.md"
    "MISMIP+ part two: retreat experiments"         => "examples/mismip_plus_retreat.md"
    "Two-dimensional flow on a bumpy bed"           => "examples/bumpy_bed.md"  
    "Melt rate parametrizations"                    => "examples/melt_parametrizations.md" 
    "Ice sheet retreat under stochastic forcing"    => "examples/stochastic_forcing.md"                       
    "Real world: West Antarctic Ice Sheet"          => "examples/WAIS.md"
#    "Two-dimensional flow with overdeepened bed" => "generated/overdeepened_bed.md" ,
#    "Melt rate parametrizations" => "generated/melt_rate_parametrizations.md" 
]


API_pages = [
    "Overview" => "API/overview.md"
    "Simulations" => "API/simulations.md"
    "Models" => "API/model.md"
    "Timestepping Parameters" => "API/timestepping_params.md"
    "Output Writing" => "API/output_writing.md"
    "Grid" => "API/grid.md"
    "Solver Parameters" => "API/solver_parameters.md"
    "Specifications" => "API/specifications.md"
    "Physical Parameters" => "API/params.md"
    "Initial Conditions" => "API/initial_conditions.md"
    "Fields" => "API/fields.md"
    "Basal Melt Rate Models" => "API/melt_rate_models.md"
]


physics_pages = [
    "Overview" => "physics/overview.md"
    "Governing Equations" => "physics/governing_equations.md"
    "Basal Melt Rate Parametrizations" => "physics/melting.md"
    #"Calving" => "physics/calving.md"
    #"Damage" => "physics/damage.md"
]

#numerics_pages = [
#    "Numerical Implementation" => "numerical_procedure/numerical_procedure.md"
#]

pages = [
    "Home" => "index.md",
    "Installation instructions" => [
        "installation_instructions.md",
        "MPI setup" => "mpi_setup.md",
        "GPU setup" => "gpu_setup.md",
    ],
    "Examples" => example_pages,
    "Running on HPC" => "running_on_hpc.md",
    "Benchmarking" => "benchmarking.md",
    "Physics" => physics_pages,
    "Numerical Implementation" => "numerical_procedure/numerical_procedure.md",
    "Model Specifications" => "model_specifications.md",
    "API" => API_pages,
   # "Simulation tips" => "simulation_tips.md",
    "MITgcm coupling" => "mitgcm_coupling.md",
    "Contributors guide" => "contributing.md",
    "Unit Testing" => "test.md",
    "Contact us" => "contact.md",
    "References" => "references.md",
]


#####
##### Build and deploy docs
#####

format = Documenter.HTML(
    collapselevel = 1,
    prettyurls = true,
    #canonical = "website_url_here",
    mathengine = MathJax3()
)


makedocs(
    sitename = "WAVI.jl",
    format = format,
    pages = pages,
    modules = [WAVI],
    doctest = ("doctest" in ARGS),
    linkcheck = ("linkcheck" in ARGS),
    checkdocs = :none,
    clean = false,
    plugins = [bib],
)


deploydocs(
    repo="github.com/WAVI-ice-sheet-model/WAVI.jl",
    devbranch="main",
    versions = nothing
)


