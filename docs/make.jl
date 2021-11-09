push!(LOAD_PATH,"../src/")
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

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUT_DIR   = joinpath(@__DIR__, "src","generated")

examples = [
    "planar_one_dimensional_flow.jl"
    "bumpy_bed.jl"
    "overdeepened_bed.jl"
    "melt_rate_parametrizations.jl"
    "west_antarctica.jl"
    #variable_slipperiness.jl
]

for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUT_DIR; flavor = Literate.DocumenterFlavor())
end


#####
#### Organize page hierarchies
#####

example_pages = [
    "One-dimensional planar flow"    => "generated/planar_one_dimensional_flow.md",
    "Flow over a bumpy bed"          => "generated/bumpy_bed.md",
    "Two-dimensional flow with overdeepened bed" => "generated/overdeepened_bed.md" ,
    "Melt rate parametrizations" => "generated/melt_rate_parametrizations.md" 
]


data_structure_pages = [
    "Overview" => "data_structure/overview.md"
    "Simulations" => "data_structure/simulations.md"
    "Models" => "data_structure/model.md"
    "Timestepping Parameters" => "data_structure/timestepping_params.md"
    "Output Writing" => "data_structure/output_writing.md"
    "Grid" => "data_structure/grid.md"
    "Solver Parameters" => "data_structure/solver_parameters.md"
    "Physical Paramters" => "data_structure/params.md"
    "Initial Conditions" => "data_structure/initial_conditions.md"
    "Fields" => "data_structure/fields.md"
    "Melt Rate Models" => "data_structure/melt_rate_models.md"
]


physics_pages = [
    "Overview" => "physics/overview.md"
    "Governing Equations" => "physics/governing_equations.md"
    "Melt Rates" => "physics/melting.md"
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
    "WAVI Setup" => data_structure_pages,
    "Simulation tips" => "simulation_tips.md",
    "MITgcm coupling" => "mitgcm_coupling.md",
    "Contributor's guide" => "contributing.md",
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
#makedocs(sitename="My Documentation")

