
#Struct to hold information on h-grid, located at cell centers.
@with_kw struct HGrid{T <: Real, N <: Integer}
    nx::N
    ny::N
    x0::T = 0.0
    y0::T = 0.0
    dx::T = 1.0
    dy::T = 1.0
    xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
    yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
    mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny); @assert mask == clip(mask)
    n::N = count(mask); @assert n == count(mask)
    crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
    samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
    spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
    b::Array{T,2} = params.bed_elevation; @assert size(b)==(nx,ny)
    h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
    s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
    dhdt::Array{T,2} = zeros(nx,ny); @assert size(dhdt)==(nx,ny)
    accumulation::Array{T,2} = zeros(nx,ny); @assert size(accumulation)==(nx,ny)
    basal_melt::Array{T,2} = zeros(nx,ny); @assert size(basal_melt)==(nx,ny)
    haf::Array{T,2} = zeros(nx,ny); @assert size(haf)==(nx,ny)
    grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
    dsdh::Array{T,2} = ones(nx,ny); @assert size(dsdh)==(nx,ny)
    shelf_strain_rate::Array{T,2} = zeros(nx,ny); @assert size(shelf_strain_rate)==(nx,ny)
    av_speed::Array{T,2} = zeros(nx,ny); @assert size(av_speed)==(nx,ny)
    bed_speed::Array{T,2} = zeros(nx,ny); @assert size(bed_speed)==(nx,ny)
    weertman_c::Array{T,2} = zeros(nx,ny); @assert size(weertman_c)==(nx,ny)
    β::Array{T,2} = zeros(nx,ny); @assert size(β)==(nx,ny)
    βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
    τbed::Array{T,2} = zeros(nx,ny); @assert size(τbed)==(nx,ny)
    ηav::Array{T,2}; @assert size(ηav)==(nx,ny)
    quad_f2::Array{T,2} = h./(3*ηav); @assert size(quad_f2)==(nx,ny)
    dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
    dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
    end
    
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

#Struct to hold information on v-grid, located at grid-North and grid-South cell faces.
@with_kw struct VGrid{T <: Real, N <: Integer}
nx::N
ny::N
x0::T = 0.0
y0::T = 0.0
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+(j-1.0)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
cent::KronType{T,N} = c(ny-1) ⊗ spI(nx)
∂x::KronType{T,N} = χ(ny-2) ⊗ ∂1d(nx-1,dx)
∂y::KronType{T,N} = ∂1d(ny-1,dy) ⊗ spI(nx)
levels::N
dwt::KronType{T,N} = wavelet_matrix(ny,levels,"forward" ) ⊗ wavelet_matrix(nx,levels,"forward")
s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(-βeff[:])*crop)
v::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
end

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