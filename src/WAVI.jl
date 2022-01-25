module WAVI

#Useful packages
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
      IterativeSolvers, Interpolations, BenchmarkTools, Reexport,
      NetCDF, JLD2, Setfield, MAT, ImageFiltering

#Import functions so they can be modified in this module.
import LinearAlgebra: ldiv!
import SparseArrays: spdiagm, spdiagm_internal, dimlub
import Setfield: @set

#This module will export these functions and types, allowing basic use of the model.
export
    #Structures
    Model, Params, TimesteppingParams, Grid, SolverParams, InitialConditions, OutputParams, Simulation,

    #Simulation controls
    update_state!, timestep!, run_simulation!,

    #Melt rates
    PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, UniformMeltFloatOnly

    #Post-processing controls
    volume_above_floatation, height_above_floatation

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2
@reexport using Setfield

#Abstract types
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractMeltRate end
abstract type AbstractModel{T <: Real, N <: Integer, M <: AbstractMeltRate} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end
#abstract type AbstractSimulation{T,N,R,A,W} end



#Type alias, just for abreviation
const KronType{T,N} = LinearMaps.KroneckerMap{T,Tuple{LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}},
                        LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}}}} where {T <: Real, N <: Integer}

#Concrete types

##################################################################################
#include all of the code
include("OutputParams/OutputParams.jl")
include("Grid.jl")
include("Params.jl")
include("SolverParams.jl")
include("TimesteppingParams.jl")
include("Clock.jl")
include("InitialConditions.jl")
include("Wavelets/Wavelets.jl")
include("Fields/Fields.jl")
include("Models/Model.jl")
include("MeltRate/MeltRate.jl")
include("Simulations/Simulation.jl")
include("utilities.jl")


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


