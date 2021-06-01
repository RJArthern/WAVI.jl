#Struct to hold information on c-grid, located at cell corners.
@with_kw struct CGrid{T <: Real, N <: Integer}
nx::N
ny::N
x0::T = 0.0
y0::T = 0.0
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+i*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+j*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
cent::KronType{T,N} = sparse(c(ny)') ⊗ sparse(c(nx)')
end
