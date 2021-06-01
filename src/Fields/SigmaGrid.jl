#Struct to hold information on 3d-grid, (extends h-grid to multiple sigma levels).
@with_kw struct SigmaGrid{T <: Real, N <: Integer}
nx::N
ny::N
nσ::N
σ::Vector{T} = collect(range(0.0,length=nσ,stop=1.0)); @assert length(σ) == nσ
ζ::Vector{T} = one(eltype(σ)) .- σ ; @assert length(ζ) == nσ
quadrature_weights::Vector{T} = [0.5;ones(nσ-2);0.5]/(nσ-1); @assert length(quadrature_weights) == nσ
η::Array{T,3} = fill(params.starting_viscosity,nx,ny,nσ); @assert size(η)==(nx,ny,nσ)
θ::Array{T,3} = fill(params.starting_temperature,nx,ny,nσ); @assert size(θ)==(nx,ny,nσ)
Φ::Array{T,3} = fill(params.starting_damage,nx,ny,nσ); @assert size(Φ)==(nx,ny,nσ)
glen_b::Array{T,3} = glen_b.(θ,Φ); @assert size(glen_b)==(nx,ny,nσ)
end