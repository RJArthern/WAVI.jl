
#Struct to hold information on wavelet-grid (u-component).
@with_kw struct UWavelets{T <: Real, N <: Integer}
nx::N
ny::N
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::Base.RefValue{N} = Ref(count(mask)); @assert n[] == count(mask)
crop::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]))
samp::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(crop[][mask[:],:]); @assert samp[] == crop[][mask[:],:]
spread::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
levels::N
idwt::KronType{T,N} = wavelet_matrix(ny,levels,"reverse" ) ⊗ wavelet_matrix(nx,levels,"reverse")
wavelets::Array{T,2} = zeros(nx,ny); @assert size(wavelets)==(nx,ny)
end

#Struct to hold information on wavelet-grid (v-component).
@with_kw struct VWavelets{T <: Real, N <: Integer}
nx::N
ny::N
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::Base.RefValue{N} = Ref(count(mask)); @assert n[] == count(mask)
crop::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]))
samp::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(crop[][mask[:],:]); @assert samp[] == crop[][mask[:],:]
spread::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
levels::N
idwt::KronType{T,N} = wavelet_matrix(ny,levels,"reverse" ) ⊗ wavelet_matrix(nx,levels,"reverse")
wavelets::Array{T,2} = zeros(nx,ny); @assert size(wavelets)==(nx,ny)
end



#Struct to hold information about wavelet-based multigrid preconditioner.
@with_kw struct Preconditioner{T <: Real, N <: Integer} <: AbstractPreconditioner{T,N}
op::LinearMap{T}
op_diag::Vector{T} = diag(sparse(op))
nsmooth::N = 5
sweep::Vector{N}
sweep_order::Vector{N} = unique(sweep)
smoother_omega::T = 1.0
restrict::LinearMap{T}
prolong::LinearMap{T}
op_coarse::LinearMap{T} = restrict*op*prolong
correction_coarse::Vector{T} = zeros(T,size(op_coarse,2))
tol_coarse::T = 1e-7
maxiter_coarse::N = 1000
end


include("preconditioners.jl")
include("update_wavelets.jl")