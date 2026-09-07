module Wavelets

using LinearAlgebra
using LinearMaps
using Parameters
using SparseArrays

using WAVI: AbstractModel, MapOrMatrix, AbstractPreconditioner
import WAVI.Utilities: fill_index_map!, stencil_scratch!, haar_dwt!, MultigridScratch
using WAVI.Stencils

include("UWavelets.jl")
include("VWavelets.jl")
include("preconditioners.jl")
include("update_preconditioners.jl")
include("update_wavelets.jl")

end