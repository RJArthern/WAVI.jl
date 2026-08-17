module Wavelets

using InplaceOps
using LinearAlgebra
using LinearMaps
using Parameters
using SparseArrays

using WAVI: AbstractModel, MapOrMatrix, AbstractPreconditioner
using WAVI.KroneckerProducts

include("UWavelets.jl")
include("VWavelets.jl")
include("preconditioners.jl")
include("update_preconditioners.jl")
include("update_wavelets.jl")

end