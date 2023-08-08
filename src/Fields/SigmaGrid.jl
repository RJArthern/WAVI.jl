"""
SigmaGrid(;
        nxs,
        nys,
        nσs,
        σ = collect(range(0.0,length=nσs,stop=1.0)),
        ζ = one(eltype(σ)) .- σ 
        quadrature_weights = [0.5;ones(nσs-2);0.5]/(nσs-1)
        η,
        θ,
        Φ,
        glen_b
        )

Construct a WAVI.jl SigmaGrid with size (nxs,nys,nσs)
SigmaGrid stores fields that are defined on the problem's three dimensional Sigma grid. 


Keyword arguments
=================
- 'nxs': (required) Number of grid cells in x-direction in SigmaGrid (should be same as grid.nx)
        Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
- 'nys': (required) Number of grid cells in y-direction in SigmaGrid (should be same as grid.ny )
- 'nσs': (required) Number of grid cells in the vertical
- 'σ' : (required) Vertical levels, read in from grid
- 'ζ' : Reverse vertical levels
- 'quadrature_weights' : (required) weights associated with sigma levels, used in computation of integrals over thickness
- 'η' : (required) three dimensional viscosity field
- 'θ' : (required) three dimensional temperature field
- 'Φ' : (required) three dimensional damage field
- 'glen_b': (required) three dimensional field of glen_b values in viscosity calcluations
"""

@with_kw struct SigmaGrid{T <: Real, N <: Integer}
nxs :: N
nys :: N
nσs :: N
σ :: Vector{T} 
ζ :: Vector{T} = one(eltype(σ)) .- σ ; @assert length(ζ) == nσs
quadrature_weights :: Vector{T} = 0.5*[ σ[2] .- σ[1] ; σ[3:end] .- σ[1:end-2] ; σ[end] .- σ[end-1] ] ; @assert length(quadrature_weights) == nσs
η :: Array{T,3}; @assert size(η)==(nxs,nys,nσs)
θ :: Array{T,3}; @assert size(θ)==(nxs,nys,nσs)
Φ :: Array{T,3}; @assert size(Φ)==(nxs,nys,nσs)
glen_b :: Array{T,3} = glen_b.(θ,Φ); @assert size(glen_b)==(nxs,nys,nσs)
end
