push!(LOAD_PATH,"../src/")
import Pkg; 
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.add(name="Documenter", version="0.27.25")
Pkg.instantiate()

using Documenter
using DocumenterCitations
using Literate 
using Plots
using WAVI

bib_filepath = joinpath(dirname(@__FILE__), "wavi.bib")
bib = CitationBibliography(bib_filepath)


ENV["GKSwstype"] = "100"

#####
##### Generate examples
#####

#const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
#const OUT_DIR   = joinpath(@__DIR__, "src","generated")

#examples = [# "planar_1D_flow.jl"
#    "bumpy_bed.jl"
#    "overdeepened_bed.jl"
#    "melt_rate_parametrizations.jl"
#    "west_antarctica.jl"
    #variable_slipperiness.jl
#]

#for example in examples
#    example_filepath = joinpath(EXAMPLES_DIR, example)
#    Literate.markdown(example_filepath, OUT_DIR; flavor = Literate.DocumenterFlavor())
#end


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
    "Installation instructions" => "installation_instructions.md",
    "Examples" => example_pages,
    "Physics" => physics_pages,
    "Numerical Implementation" => "numerical_procedure/numerical_procedure.md",
    "API" => API_pages,
   # "Simulation tips" => "simulation_tips.md",
    "MITgcm coupling" => "mitgcm_coupling.md",
    "Contributors guide" => "contributing.md",
    "Contact us" => "contact.md",
    "References" => "references.md",
]


#####
##### Build and deploy docs
#####

format = Documenter.HTML(
    collapselevel = 1,
       prettyurls = get(ENV, "CI", nothing) == "true",
        #canonical = "website_url_here",
       mathengine = MathJax3()
)


makedocs(bib,
    sitename = "WAVI.jl",
    format = format,
    pages = pages,
    modules = [WAVI],
    doctest = false,
    strict = false,
    clean = false,
    checkdocs = :none)

deploydocs(
    repo="github.com/RJArthern/WAVI.jl",
    devbranch="main",
    versions = nothing
)


