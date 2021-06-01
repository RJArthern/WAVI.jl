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
u::Array{T,2} = zeros(nx,ny); @assert size(u)==(nx,ny)
v::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
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
