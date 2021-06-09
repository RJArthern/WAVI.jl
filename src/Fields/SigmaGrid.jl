"""
SigmaGrid(;
        Nx,
        Ny,
        Nσ,
        σ = collect(range(0.0,length=Nσ,stop=1.0)),
        ζ = one(eltype(σ)) .- σ 
        quadrature_weights = [0.5;ones(Nσ-2);0.5]/(Nσ-1)
        η,
        θ,
        Φ,
        glen_b
        )

Construct a WAVI.jl SigmaGrid with size (Nx,Ny,Nσ)
SigmaGrid stores fields that are defined on the problem's three dimensional Sigma grid. 


Keyword arguments
=================
- 'Nx': (required) Number of grid cells in x-direction in SigmaGrid (should be same as grid.nx)
        Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
- 'Ny': (required) Number of grid cells in y-direction in SigmaGrid (should be same as grid.ny )
- 'Nσ': (required) Number of grid cells in the vertical
- 'σ' : Vertical levels
- 'ζ' : Reverse vertical levels
- 'quadrature_weights' : (required) weights associated with sigma levels, used in computation of integrals over thickness
- 'η' : (required) three dimensional viscosity field
- 'θ' : (required) three dimensional temperature field
- 'Φ' : (required) three dimensional damage field
- 'glen_b': (required) three dimensional field of glen_b values in viscosity calcluations
"""

@with_kw struct SigmaGrid{T <: Real, N <: Integer}
Nx :: N
Ny :: N
Nσ :: N
σ :: Vector{T} = collect(range(0.0,length=Nσ,stop=1.0)); @assert length(σ) == Nσ
ζ :: Vector{T} = one(eltype(σ)) .- σ ; @assert length(ζ) == Nσ
quadrature_weights :: Vector{T} = [0.5;ones(Nσ-2);0.5]/(Nσ-1); @assert length(quadrature_weights) == Nσ
η :: Array{T,3}; @assert size(η)==(Nx,Ny,Nσ)
θ :: Array{T,3}; @assert size(θ)==(Nx,Ny,Nσ)
Φ :: Array{T,3}; @assert size(Φ)==(Nx,Ny,Nσ)
glen_b :: Array{T,3} = glen_b.(θ,Φ); @assert size(glen_b)==(Nx,Ny,Nσ)
end
