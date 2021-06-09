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
    
    
    
    