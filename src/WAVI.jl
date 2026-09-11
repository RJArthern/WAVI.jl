module WAVI

# # TODO Needs tidy up. Refactoring has moved (or will move) many of these into the submodules where they are used
using Reexport

# #Import functions so they can be modified in this module.
# Left commented for now because I don't think we use this here anymore.
# import Base: *, size, eltype
# import LinearAlgebra: ldiv!,mul!
import Setfield: @set

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2
@reexport using Setfield

# Abstract types

#Numerical strategy
abstract type AbstractSpec end

#Core functionality
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractField{T <: Real, N <: Integer} end
abstract type AbstractModel{T <: Real,
                            N <: Integer,
                            S <: AbstractSpec,
                            F <: AbstractField,
                            G <: AbstractGrid} end

#Momentum and thickness solves and simulation control                    
abstract type AbstractPreconditioner{T <: Real,N <: Integer} end
abstract type AbstractSimulation end
abstract type AbstractClimateForcing end
abstract type AbstractInversionSimulation end


#Modular physics
abstract type AbstractMeltRate end
abstract type AbstractSurfaceMassBalance end
abstract type AbstractFracture end
abstract type AbstractSlidingLaw end
abstract type AbstractBasalHydrology end
abstract type AbstractThermoDynamics end

using LinearMaps
#Type alias, just for abreviation
const MapOrMatrix{T} = Union{LinearMap{T}, AbstractMatrix{T}}

##################################################################################
#include all of the code
include("Deferred.jl")
include("Time.jl")
include("Grids.jl")
include("Parameters.jl")
include("ClimateForcing/ClimateForcing.jl")
include("KroneckerProducts.jl")
include("Utilities.jl")
include("Wavelets/Wavelets.jl")
include("Fields/Fields.jl")
include("MeltRates/MeltRates.jl")
include("SurfaceMassBalance/SurfaceMassBalance.jl")
include("Advection/Advection.jl")
include("Fracture/Fracture.jl")
include("SlidingLaw/SlidingLaw.jl")
include("BasalHydrology/BasalHydrology.jl")
include("ThermoDynamics/ThermoDynamics.jl")
include("Processes/Processes.jl")
include("Models/Models.jl")
include("Outputs/Outputs.jl")
include("Simulations/Simulation.jl")
include("Specs/Specs.jl")
include("Inversion/Inversion.jl")


export AbstractField, AbstractGrid, AbstractMeltRate, AbstractSurfaceMassBalance,
  AbstractFracture, AbstractSlidingLaw , AbstractBasalHydrology,
   AbstractThermoDynamics, AbstractModel, AbstractPreconditioner,
   AbstractSpec,AbstractClimateForcing

using .Deferred
export Collector, clear!, collect!, register_field!

using .Time
export Clock, compute_iterations_and_end_time

using .Parameters
export Params, SolverParams, TimesteppingParams

using .KroneckerProducts
export KronType, KroneckerProduct

using .Grids
export Grid

using .Utilities
export volume_above_floatation, height_above_floatation, get_c_mask, clip

using .Wavelets
export Preconditioner 

using .Fields
export GridField, InitialConditions

using .MeltRates
export PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, 
QuadraticForcedMeltRate, MeltRateExponentVariation, MeltRateExponentVariationBasins, UniformMeltUnderShelves, 
UniformMeltUnderShelvesBasins, ISMIP7MeltRate

using .SurfaceMassBalance
export AccumulationFromParams, BinfileAccumulationRate, UniformAccumulationRate, ISMIP7SMB
using .Processes
export update_state!, update_velocities!

using .Models
export AbstractModel, Model, update_state!

using .Outputs
export OutputParams,
    fetch_output, 
    get_spatiotemporal_var_atts, get_spatial_dimensions, get_times, get_output_as_dict, 
    make_ncfile, make_ncfile_from_filenames, 
    write_output, write_outputs, zip_output

using .Simulations
export Simulation, run_simulation!, timestep!, 
    # TODO: these probably should be in clock, processes and outputs?
    update_clock!, update_thickness!, write_vel

using .Specs
export BasicSpec, ThreadedSpec, MPISpec

using .Fracture
export ConstantDamage, DruckerPragerPhaseField, ISMIP7Hydrofracture

using .SlidingLaw
export WeertmanSlidingLaw, WeertmanRelaxationSlidingLaw, CoulombSlidingLaw, BuddSlidingLaw, TsaiSlidingLaw, TsaiBuddSlidingLaw, SchoofSlidingLaw, ZoetIversonSlidingLaw

using .BasalHydrology
export NoHydrology, ConstantBasalWaterThickness, LeakyBucket, SheetOnlyGlaDS

using .ThermoDynamics
export NoThermoDynamics, QuadraticTemperatureApproximation, QuadraticTemperatureApproximationIcebergTest

using .Inversions
export run_inversion_simulation!
export Inversion, InversionParams, JKVsteppingParams, DataFields, 
       InversionSimulation, InversionOutput, smooth_veloc_data

using .ClimateForcing 
export ISMIP7_ANOMALY,ISMIP7_CONTROL,ISMIP7_OCX

end


