push!(LOAD_PATH,"../src/")
#import Pkg; Pkg.add("Documenter")
#import Pkg; Pkg.add("DocumenterCitations")
#import Pkg; Pkg.add("Literate")
#import Pkg; Pkg.add("Plots")

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


data_structure_pages = [
    "Overview" => "data_structure/overview.md"
    "Simulations" => "data_structure/simulations.md"
    "Models" => "data_structure/model.md"
    "Timestepping Parameters" => "data_structure/timestepping_params.md"
    "Output Writing" => "data_structure/output_writing.md"
    "Grid" => "data_structure/grid.md"
    "Solver Parameters" => "data_structure/solver_parameters.md"
    "Physical Parameters" => "data_structure/params.md"
    "Initial Conditions" => "data_structure/initial_conditions.md"
    "Fields" => "data_structure/fields.md"
    "Melt Rate Models" => "data_structure/melt_rate_models.md"
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
    "API" => data_structure_pages,
   # "Simulation tips" => "simulation_tips.md",
    "MITgcm coupling" => "mitgcm_coupling.md",
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
   authors = "Alexander Bradley",
    format = format,
     pages = pages,
   modules = [WAVI],
   doctest = false,
    strict = false,
     clean = false,
 checkdocs = :none # Should fix our docstring so we can use checkdocs=:exports with strict=true.
)

deploydocs(;
    repo="github.com/RJArthern/WAVI.jl",
    devbranch="build-docs",
    versions = nothing
)


