module Models

using Parameters
using Setfield

using WAVI: AbstractField, 
            AbstractGrid, 
            AbstractMeltRate,
            AbstractSurfaceMassBalance, 
            AbstractFracture,
            AbstractSlidingLaw,
            AbstractBasalHydrology,
            AbstractThermoDynamics,
            AbstractModel, 
            AbstractSpec
using WAVI.Architectures: on_architecture, CPU
import WAVI.Architectures: architecture
using WAVI.Fields
using WAVI.Grids
using WAVI.MeltRates
using WAVI.SurfaceMassBalance
using WAVI.Fracture
using WAVI.SlidingLaw
using WAVI.BasalHydrology
using WAVI.ThermoDynamics
using WAVI.Parameters

export Model, update_state!, restore_pickup_architecture!

"""
    BasicSpec()

Default specification: one process, arrays on the CPU.

If no `spec` is passed to Model, a `BasicSpec` is used. Stencil
kernels still use Julia threads when you launch with `julia -t N`. This
is different to `ThreadedSpec`, which splits the domain into overlapping
Schwarz subdomains which I think should be removed (or, not advertised
in the docs) in the future.

TODO: Remove ThreadedSpec in a future optimisation.
"""
struct BasicSpec <: AbstractSpec 
    function BasicSpec()
        @debug "Not implementing any parallel computations, running with BasicSpec"
        return new()
    end
end

struct Model{T,N,S,F,G,
             M<:AbstractMeltRate,
             SMB<:AbstractSurfaceMassBalance,
             FR<:AbstractFracture,
             SL<:AbstractSlidingLaw,
             BH<:AbstractBasalHydrology,
             TD<:AbstractThermoDynamics
             } <: AbstractModel{T,N,S,F,G}
    grid    ::  G
    fields  ::  F
    params  ::  Params
    solver_params :: SolverParams
    spec   ::  S
    shelf_melt_rate :: M
    surface_mass_balance :: SMB
    fracture :: FR
    sliding_law :: SL
    basal_hydrology :: BH
    thermo_dynamics :: TD
    verbose :: Bool
end

"""
    Model()

    Construct a WAVI.jl model object.

"""
function Model(grid::G, 
               bed_elevation::Union{Integer, Function, AbstractArray}, 
               spec::S;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               shelf_melt_rate::M = UniformMeltRate(),
               surface_mass_balance::SMB = AccumulationFromParams(),
               fracture::FR = ConstantDamage(),
               sliding_law::SL = WeertmanSlidingLaw(),
               basal_hydrology::BH = NoHydrology(),
               thermo_dynamics::TD = NoThermoDynamics(),
               verbose = true
               ) where {G<:AbstractGrid, 
                        S<:AbstractSpec, 
                        M<:AbstractMeltRate,
                        SMB<:AbstractSurfaceMassBalance,
                        FR<:AbstractFracture,
                        SL<:AbstractSlidingLaw,
                        BH<:AbstractBasalHydrology,
                        TD<:AbstractThermoDynamics}

    spec_label = hasproperty(spec, :ngridsx) ?
        "$(nameof(typeof(spec))) $(spec.ngridsx)x$(spec.ngridsy)" :
        string(nameof(typeof(spec)))
    @debug "$(spec_label) grid $(grid.nx)x$(grid.ny)"

    # FIXME: this all smells, hacking for threading
    bed_array = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    
    #expand scalar paramaters onto grid
    params = reconstruct_on_grid(params,grid)

    #Replace all NaN entries with defaults from params on correct grid              
    initial_conditions = reconstruct_on_grid(initial_conditions, params, grid)

    #expand any spatial parameters needed onto correct grid
    shelf_melt_rate = reconstruct_on_grid(shelf_melt_rate,grid)
    surface_mass_balance = reconstruct_on_grid(surface_mass_balance,grid)
    fracture = reconstruct_on_grid(fracture,grid)
    sliding_law = reconstruct_on_grid(sliding_law,grid)
    basal_hydrology = reconstruct_on_grid(basal_hydrology,grid)
    thermo_dynamics = reconstruct_on_grid(thermo_dynamics,grid)

    
    # TODO: grids are heavily reliant on the use of keyword arguments which do not allow specializations / multiple dispatch to work effectively

    # TODO: the passthrough of arguments like this is smelly - Configuration should be a type
    fields = GridField(grid, bed_array; initial_conditions, params, solver_params)
    arch = architecture(spec)
    if !(arch isa CPU)
        fields = on_architecture(arch, fields)
    end

    model = Model(
               grid, 
               fields, 
               params, 
               solver_params, 
               spec, 
               shelf_melt_rate,
               surface_mass_balance,
               fracture,
               sliding_law,
               basal_hydrology,
               thermo_dynamics,
               verbose)
    @debug "$(nameof(typeof(spec))): model ready"
    return model
end

Model(grid, bed_elev; kw...) = Model(grid, bed_elev, BasicSpec(); kw...)
Model(; grid, bed_elevation, spec = BasicSpec(), kw...) = Model(grid, bed_elevation, spec; kw...)

architecture(model::AbstractModel) = architecture(model.spec)

"""
    restore_pickup_architecture!(model)

Put dense fields on the spec's architecture and drop stencil scratch.

Call this after loading a checkpoint, before the next velocity solve.
JLD2 may restore mixed host and device arrays, or a host scratch buffer,
while the spec says the run is on the GPU. The file format is left
unchanged, so later checkpoint layouts can keep calling this hook.
Do not use this on MPI `spec.global_fields`; those stay on the host.
"""
function restore_pickup_architecture!(model::AbstractModel)
    model.fields.stencil_scratch[] = nothing
    arch = architecture(model)
    arch isa CPU && return model
    fields = on_architecture(arch, model.fields)
    return @set model.fields = fields
end

# This is to enable use of Setfield, which derives a parameter setup from the fields of an existing structure via JuliaObjects
# FIXME: this wasn't required in the original WAVI codebase. 
#   Model needs to be in a position to have type analysis on it's properties to recreate the instance of it by ConstructionBase, which requires some refactoring
#   Ref: https://juliaobjects.github.io/ConstructionBase.jl/dev/#type-tips
Model(g::G, f::F, p::P, sp::SP, s::S, m::M, smb::SMB, fr::FR, sl::SL, bh::BH, td::TD, vb::Bool
    ) where {
    G<:AbstractGrid,
    F<:AbstractField,
    P<:Params,
    SP<:SolverParams,
    S<:AbstractSpec,
    M<:AbstractMeltRate,
    SMB<:AbstractSurfaceMassBalance,
    FR<:AbstractFracture,
    SL<:AbstractSlidingLaw,
    BH<:AbstractBasalHydrology,
    TD<:AbstractThermoDynamics
    } = 
    Model{Float64,Int64,S,F,G,M,SMB,FR,SL,BH,TD}(
                            g, f, p, sp, s, m, smb, fr, sl, bh, td, vb)

##
# Global domain alterations
#
Base.propertynames(model::Model{T,N,S,F,G,M,SMB,FR,SL,BH,TD}, private::Bool) where {T,N,S,F,G,M,SMB,FR,SL,BH,TD} = (fieldnames(typeof(model))..., :global_fields, :global_grid)

# FIXME: this is a sign of a frustration in WAVIs structural layout - too many deep nested structures accessed through high level passing
#  which inhibits multiple dispatch
function Base.getproperty(model::Model{T,N,S,F,G,M,SMB,FR,SL,BH,TD}, s::Symbol) where {T,N,S,F,G,M,SMB,FR,SL,BH,TD}
    if s == :global_fields
        return getfield(model, :fields)
    elseif s == :global_grid
        return getfield(model, :grid)
    end
    getfield(model, s)
end


include("utils.jl")

function Base.show(io::IO, m::Model)
    return print(io, "Characteristics of ", summary(m))
end

end
