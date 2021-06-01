module WAVI

#Useful packages
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
      IterativeSolvers, Interpolations, BenchmarkTools, Reexport,
      NetCDF, JLD2, HDF5, Setfield, MAT

#Import functions so they can be modified in this module.
import LinearAlgebra: ldiv!
import SparseArrays: spdiagm, spdiagm_internal, dimlub
import Setfield: @set

#This module will export these functions and types, allowing basic use of the model.
export update_state!, timestep!, Model, Params, TimesteppingParams, 
Grid, SolverParams, InitialConditions, OutputParams, Simulation, run_simulation!

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2
@reexport using Setfield

#Abstract types
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractModel{T <: Real, N <: Integer} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end
abstract type AbstractSimulation{T <: Real, N <: Integer, R <: Real} end
abstract type AbstractMeltRateModel{PC <: Bool, M} end


#Type alias, just for abreviation
const KronType{T,N} = LinearMaps.KroneckerMap{T,Tuple{LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}},
                        LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}}}} where {T <: Real, N <: Integer}

#Concrete types

##################################################################################
#type definitions
include("./OutputParams/OutputParams.jl")
include("./Grid.jl")
include("./Params.jl")
include("./Clock.jl")
include("./InitialConditions.jl")
include("./Fields/HGrid.jl")
include("./Fields/CGrid.jl")
include("./Fields/UGrid.jl")
include("./Fields/VGrid.jl")
include("./Fields/SigmaGrid.jl")
include("./WaveletGrids.jl")
include("./Fields.jl")
include("./Model.jl")
include("./MeltRateModels/MeltRateModel.jl")
include("./Simulation.jl")




#Functions
include("./OutputParams/output_writing.jl")
include("./OutputParams/zipping_output.jl")
include("preconditioners.jl")
include("update_state.jl")
include("update_velocities.jl")
include("utilities.jl")
include("wavelets.jl")
include("run_simulation.jl")


"""
    spdiagm(m::Integer, n::Integer, kv::Pair{<:Integer,<:AbstractVector}...)

Method from Julia V1.4.2 to create non-square sparse matrix from diagonals. Included for backward compatibility.
"""
spdiagm(m::Integer, n::Integer, kv::Pair{<:Integer,<:AbstractVector}...) = _spdiagm((Int(m),Int(n)), kv...)
function _spdiagm(size, kv::Pair{<:Integer,<:AbstractVector}...)
    I, J, V = spdiagm_internal(kv...)
    mmax, nmax = dimlub(I), dimlub(J)
    mnmax = max(mmax, nmax)
    m, n = something(size, (mnmax,mnmax))
    (m ≥ mmax && n ≥ nmax) || throw(DimensionMismatch("invalid size=$size"))
    return sparse(I, J, V, m, n)
end
end


