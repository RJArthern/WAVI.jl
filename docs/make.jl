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
    "MISMIP+"   => "examples/mismip_plus.md"
]

data_structure_pages = [
    "Overview" => "data_structure/overview.md"
    "Simulations" => "data_structure/simulations.md"
    "Models" => "data_structure/model.md"
]

physics_pages = [
    "Governing equations" => "physics/governing_equations.md"
]

parametrization_pages = [
    "Overview" => "parametrizations/overview.md"
]

pages = [
    "Home" => "index.md",
    "Installation instructions" => "installation_instructions.md",
    "Data structure" => data_structure_pages,
    "Examples" => example_pages,
    "Physics" => physics_pages,
    "Parametrizations" => parametrization_pages,
    "Simulation tips" => "simulation_tips.md",
    "Contributor's guide" => "contributing.md",
    "Gallery" => "gallery.md",
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

