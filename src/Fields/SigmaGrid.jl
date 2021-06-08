#Struct to hold information on 3d-grid, (extends h-grid to multiple sigma levels).
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

#struct SigmaGrid{T <: Real, N <: Integer}
#                    Nx :: N
#                    Ny :: N
#                    Nσ :: N
#                     σ :: Vector{T} 
#                     ζ :: Vector{T} 
#  quadrature_weights :: Vector{T} 
#                    η :: Array{T,3} 
#                     θ :: Array{T,3} 
#                     Φ :: Array{T,3} 
#                glen_b :: Array{T,3}   #
#end

#function SigmaGrid(;
#                    Nx,
#                    Ny,
#                    Nσ,
#                    quadrature_weights = [0.5;ones(Nσ-2);0.5]/(Nσ-1),
#                    η,
#                    θ,
#                    Φ,
#                    glen_b)

#σ = collect(range(0.0,length=Nσ,stop=1.0)); @assert length(σ) == Nσ
#ζ = one(eltype(σ)) .- σ ; @assert length(ζ) == Nσ

#return SigmaGrid(
#                Nx,
#                Ny,
#                Nσ,
#                σ,
#                ζ,
#                quadrature_weights,
#                η, 
#                θ,
#                Φ,
#                glen_b)
#end