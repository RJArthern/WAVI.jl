module Specs

using WAVI: AbstractSpec

export AbstractDecompSpec

abstract type AbstractDecompSpec <: AbstractSpec end

include("utils.jl")
include("basic.jl")
include("threaded.jl")
include("gpu.jl")
include("mpi.jl")

end