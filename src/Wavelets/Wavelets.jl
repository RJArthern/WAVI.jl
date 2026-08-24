module Wavelets

using Adapt
using LinearAlgebra
using LinearMaps
using Parameters
using SparseArrays

using WAVI.Architectures: adapt_structure_fields
using WAVI: AbstractModel, MapOrMatrix, AbstractPreconditioner
import WAVI.Utilities: fill_index_map!, stencil_scratch!, haar_dwt!, MultigridScratch, zeros_like, copy_like
using WAVI.Stencils

include("UWavelets.jl")
include("VWavelets.jl")
include("preconditioners.jl")
include("update_preconditioners.jl")
include("update_wavelets.jl")

Adapt.adapt_structure(to, w::Union{UWavelets, VWavelets}) = adapt_structure_fields(to, w)

end