push!(LOAD_PATH,"../src/")
using Documenter
using DocumenterCitations

using  WAVI

bib_filepath = joinpath(dirname(@__FILE__), "wavi.bib")
bib = CitationBibliography(bib_filepath)

#####
##### Generate examples (uncomment for example generation -- see https://github.com/CliMA/Oceananigans.jl/blob/master/docs/make.jl)
#####

#const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
#const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")
#examples = [
#    "MISMIP_PLUS.jl"
#]

#for example in examples
#    example_filepath = joinpath(EXAMPLES_DIR, example)
#    Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
#end

#####
#### Organize page hierarchies
#####

#example_pages = [
#    "One-dimensional diffusion"          => "generated/one_dimensional_diffusion.md",
#]
example_pages = [
    "MISMIP"   => "examples/mismip.md"
    "MISMIP3D"   => "examples/mismip_3d.md"
    "MISMIP+"   => "examples/mismip_plus.md"
    "MISI"   => "examples/MISI.md"
    "Greenland"   => "examples/Greenland.md"
    "WAIS"   => "examples/WAIS.md"
    "Antarctica" => "examples/Antarctica.md"

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

numerics_pages = [
    "Numerical Grid" => "numerics/numerical_grid.md"
    "Discretization" => "numerics/discretization.md"
    "Velocity Solve" => "numerics/velocity_solve.md"
    "Timestepping" => "numerics/timestepping.md"
]

physics_pages = [
    "Overview and Background" => "physics/overview.md"
    "Governing Equations" => "physics/governing_equations.md"
    "Basal Stress" => "physics/basal_stress.md"
    "Melt Rates" => "physics/melting.md"
    #"Calving" => "physics/calving.md"
    "Damage" => "physics/damage.md"
]


pages = [
    "Home" => "index.md",
    "Installation instructions" => "installation_instructions.md",
    "Examples" => example_pages,
    "Physics" => physics_pages,
    "Numerical Implementation" => numerics_pages,
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

