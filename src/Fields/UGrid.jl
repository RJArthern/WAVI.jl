#Struct to hold information on u-grid, located at grid-East and grid-West cell faces.
@with_kw struct UGrid{T <: Real, N <: Integer}
    nx::N
    ny::N
    x0::T = 0.0
    y0::T = 0.0
    dx::T = 1.0
    dy::T = 1.0
    xx::Array{T,2}=[x0+(i-1.0)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
    yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
    mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
    n::N = count(mask); @assert n == count(mask)
    crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
    samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
    spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
    cent::KronType{T,N} = spI(ny) ⊗ c(nx-1)
    ∂x::KronType{T,N} = spI(ny) ⊗ ∂1d(nx-1,dx)
    ∂y::KronType{T,N} = ∂1d(ny-1,dy) ⊗ χ(nx-2)
    levels::N
    dwt::KronType{T,N} = wavelet_matrix(ny,levels,"forward" ) ⊗ wavelet_matrix(nx,levels,"forward")
    s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
    h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
    grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
    βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
    dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(-βeff[:])*crop)
    u::Array{T,2} = zeros(nx,ny); @assert size(u)==(nx,ny)
    end
    