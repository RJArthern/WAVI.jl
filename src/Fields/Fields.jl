module Fields

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
export reconstruct_on_grid, reconstruct_on_subdomain
export GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid

using LinearAlgebra
using LinearMaps
using Parameters
using Setfield          # TODO: InitialConditions using this, bit of an anti-pattern?
using SparseArrays

using WAVI: AbstractField, AbstractGrid
using WAVI.Grids
using WAVI.KroneckerProducts
using WAVI.Parameters
using WAVI.Utilities
using WAVI.Wavelets


"""Trivial Kronecker product for assembly-buffer grids that never apply operators."""
_storage_only_kron() = spzeros(Float64, 1, 1) ⊗ spzeros(Float64, 1, 1)

include("UGrid.jl")
include("VGrid.jl")
include("HGrid.jl")
include("CGrid.jl")
include("SigmaGrid.jl")
include("InitialConditions.jl")
include("utils.jl")

"""
    Structure to hold all field variables in WAVI.jl
"""
struct GridField{T <: Real, N <: Integer} <: AbstractField{T, N}
    gh  :: HGrid{T,N}
    gu  :: UGrid{T,N}
    gv  :: VGrid{T,N}
    gc  :: CGrid{T,N}
    g3d :: SigmaGrid{T,N}
    wu  :: UWavelets{T,N}
    wv  :: VWavelets{T,N}
    stencil_scratch :: Ref{Any}
end

"""
    setup_fields(grid, initial_conditions, solver_params, params, bed_array)

Acts as a constructor for the fields (no explicit constructor as fields `only ever called when setting up a model)
"""

function GridField(grid::AbstractGrid, bed_array;
                   initial_conditions::InitialConditions = InitialConditions(),
                   params::Params = Params(),
                   solver_params::SolverParams = SolverParams(),
                   mpi_rank::Array{Float64, 2} = -ones(grid.nx, grid.ny),
                   assembly_buffer::Bool = false,
                   )

   
    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    h_mask = grid.h_mask 
    u_mask = get_u_mask(h_mask)
    v_mask = get_v_mask(h_mask)
    c_mask = get_c_mask(h_mask)

    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[grid.u_iszero].=false
    v_mask[grid.v_iszero].=false

    # Assembly buffers (e.g. MPISpec global_fields) only need correctly sized writable
    # arrays for gather/Bcast. Skip IC deepcopies, glen_b fill, and sparse/wavelet operators.
    if assembly_buffer
        gh=HGrid(
            nxh=grid.nx,
            nyh=grid.ny,
            mask=h_mask,
            h_isfixed = grid.h_isfixed,
            hyd_potential_isfixed = grid.hyd_potential_isfixed,
            b = bed_array,
            h = zeros(grid.nx, grid.ny),
            ηav = zeros(grid.nx, grid.ny),
            grounded_fraction = zeros(grid.nx, grid.ny),
            basal_water_thickness = zeros(grid.nx, grid.ny),
            hydraulic_potential_b = zeros(grid.nx, grid.ny),
            effective_pressure = zeros(grid.nx, grid.ny),
            basal_melt = zeros(grid.nx, grid.ny),
            θ_ave = zeros(grid.nx, grid.ny),
            preBfactor = ones(grid.nx, grid.ny),
            mpi_rank = mpi_rank,
            storage_only = true,
        )
        gu=UGrid(
            nxu=grid.nx+1,
            nyu=grid.ny,
            dx=grid.dx,
            dy=grid.dy,
            mask=u_mask,
            u_isfixed=grid.u_isfixed,
            u=zeros(grid.nx+1, grid.ny),
            levels=solver_params.levels,
            storage_only = true,
        )
        gv=VGrid(
            nxv=grid.nx,
            nyv=grid.ny+1,
            dx=grid.dx,
            dy=grid.dy,
            mask=v_mask,
            v_isfixed=grid.v_isfixed,
            v=zeros(grid.nx, grid.ny+1),
            levels=solver_params.levels,
            storage_only = true,
        )
        gc=CGrid(
            nxc=grid.nx-1,
            nyc=grid.ny-1,
            mask=c_mask,
            storage_only = true,
        )
        g3d=SigmaGrid(
            nxs=grid.nx,
            nys=grid.ny,
            nσs=grid.nσ,
            σ =grid.σ,
            η = zeros(grid.nx, grid.ny, grid.nσ),
            θ = zeros(grid.nx, grid.ny, grid.nσ),
            Φ = zeros(grid.nx, grid.ny, grid.nσ),
            strain_history = zeros(grid.nx, grid.ny, grid.nσ),
            glen_b = zeros(grid.nx, grid.ny, grid.nσ),
            quadrature_weights = grid.quadrature_weights
        )
        wu=UWavelets(nxuw=grid.nx+1,nyuw=grid.ny,levels=solver_params.levels, storage_only=true)
        wv=VWavelets(nxvw=grid.nx,nyvw=grid.ny+1,levels=solver_params.levels, storage_only=true)
        return GridField(gh,gu,gv,gc,g3d,wu,wv,Ref{Any}(nothing))
    end

    #h-grid
    #gh=HGrid(grid, params) #uncomment if using the explicit constructor method
    h =  deepcopy(initial_conditions.initial_thickness)
    grounded_fraction =  deepcopy(initial_conditions.initial_grounded_fraction)
    basal_water_thickness = deepcopy(initial_conditions.initial_basal_water_thickness)
    hydraulic_potential_b = deepcopy(initial_conditions.initial_hydraulic_potential_b)
    effective_pressure = deepcopy(initial_conditions.initial_effective_pressure)
    θ_ave = deepcopy(initial_conditions.initial_θ_ave)
    preBfactor = deepcopy(initial_conditions.initial_preBfactor)
    ηav = deepcopy(initial_conditions.initial_viscosity[:,:,1]) #set to the viscosity on the first level for now
    
    gh=HGrid(
    nxh=grid.nx,
    nyh=grid.ny,
    mask=h_mask,
    h_isfixed = grid.h_isfixed,
    hyd_potential_isfixed = grid.hyd_potential_isfixed,
    b = bed_array,
    h = h,
    ηav = ηav,
    grounded_fraction = grounded_fraction,
    basal_water_thickness = basal_water_thickness,
    hydraulic_potential_b = hydraulic_potential_b,
    effective_pressure = effective_pressure,
    θ_ave = θ_ave,
    preBfactor = preBfactor,
    mpi_rank = mpi_rank
    )

    #u-grid
    gu=UGrid(
        nxu=grid.nx+1,
        nyu=grid.ny,
        dx=grid.dx,
        dy=grid.dy,
        mask=u_mask,
        u_isfixed=grid.u_isfixed,
        u=deepcopy(initial_conditions.initial_u_veloc),
        levels=solver_params.levels
    )

    #v-grid
    gv=VGrid(
        nxv=grid.nx,
        nyv=grid.ny+1,
        dx=grid.dx,
        dy=grid.dy,
        mask=v_mask,
        v_isfixed=grid.v_isfixed,
        v=deepcopy(initial_conditions.initial_v_veloc),
        levels=solver_params.levels
    )

    #c-grid
    gc=CGrid(
        nxc=grid.nx-1,
        nyc=grid.ny-1,
        mask=c_mask
    )

    #3D-grid
    η = deepcopy(initial_conditions.initial_viscosity)
    θ = deepcopy(initial_conditions.initial_temperature)
    Φ = deepcopy(initial_conditions.initial_damage)
    strain_history = deepcopy(initial_conditions.initial_strain_history)


    g3_glen_b = zeros(size(η))
    @debug "Sigma-grid glen_b setup inputs" inputs = (
        size_g3_glen_b = size(g3_glen_b),
        size_θ = size(θ),
        size_Φ = size(Φ),
        size_glen_a_ref = size(params.glen_a_ref),
        nx = grid.nx,
        ny = grid.ny,
        nσ = grid.nσ,
        glen_a_ref = params.glen_a_ref,
    )
    for i = 1:grid.nx
        for j = 1:grid.ny
            for k = 1:grid.nσ
                g3_glen_b[i,j,k] = glen_b(θ[i,j,k],Φ[i,j,k],params.glen_a_ref[i,j], params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const)
            end
        end
    end

    g3d=SigmaGrid(
        nxs=grid.nx,
        nys=grid.ny,
        nσs=grid.nσ,
        σ =grid.σ,
        η = η,
        θ = θ,
        Φ = Φ,
        strain_history = strain_history,
        glen_b = g3_glen_b,
        quadrature_weights = grid.quadrature_weights
    )

    #Wavelet-grid, u-component.
    wu=UWavelets(nxuw=grid.nx+1,nyuw=grid.ny,levels=solver_params.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nxvw=grid.nx,nyvw=grid.ny+1,levels=solver_params.levels)
    return GridField(gh,gu,gv,gc,g3d,wu,wv,Ref{Any}(nothing))
end

end