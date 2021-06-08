#Struct to hold information on v-grid, located at grid-North and grid-South cell faces.
#@with_kw struct VGrid{T <: Real, N <: Integer}
#    nx::N
#   ny::N
#    x0::T = 0.0
#    y0::T = 0.0
#    dx::T = 1.0
#    dy::T = 1.0
#    xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
#    yy::Array{T,2}=[y0+(j-1.0)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
#    mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
#    n::N = count(mask); @assert n == count(mask)
#    crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
#    samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
#    spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
#    cent::KronType{T,N} = c(ny-1) ⊗ spI(nx)
#    ∂x::KronType{T,N} = χ(ny-2) ⊗ ∂1d(nx-1,dx)
#    ∂y::KronType{T,N} = ∂1d(ny-1,dy) ⊗ spI(nx)
#    levels::N
#    dwt::KronType{T,N} = wavelet_matrix(ny,levels,"forward" ) ⊗ wavelet_matrix(nx,levels,"forward")
#    s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
#    h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
#    grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
#    βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
#    dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(-βeff[:])*crop)
#    v::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#    end

struct VGrid{T <: Real, N <: Int}
                   Nx :: N                                     # Number of frid cells in x-direction in UGrid
                   Ny :: N                                     # Number of grid cells in y-direction in UGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
                    n :: N                                     # Total number of cells in model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
                 cent :: KronType{T,N}                         # Map from U grid to H grid 
                   ∂x :: KronType{T,N}                         # Matrix representation of differentiation wrt x 
                   ∂y :: KronType{T,N}                         # Matrix representation of differentiation wrt y
               levels :: N                                     # Number of levels in the vertical
                  dwt :: KronType{T,N}                         # Wavelet matrix product on u grid 
                    s :: Array{T,2}                            # Ice surface elevation
                    h :: Array{T,2}                            # Ice thickness
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction
                 βeff :: Array{T,2}                            # Effective β value on u grid(eqn 12 in Arthern 2015 JGeophysRes)
             dnegβeff :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (diagonal of βeff) 
                    v :: Array{T,2}
end

    
"""
    VGrid(;
            Nx,
            Ny,
            mask = trues(Nx,Ny), 
            levels,
            dx,
            dy)

Construct a WAVI.jl VGrid with size (Nx,Ny)
UGrid stores fields that are defined on the problem's U grid. 
(Co-ordinates of HGrid stored in a Grid under xxhu yyu fields)

Keyword arguments
=================
    - 'Nx': (required) Number of grid cells in x-direction in UGrid (should be same as grid.nx)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'Ny': (required) Number of grid cells in y-direction in HGrid (should be same as grid.ny + 1)
    - 'mask': Mask specifying the model domain with respect to U grid
    - levels: (required) Number of vertical levels 
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function VGrid(;
        Nx,
        Ny,
        mask = trues(Nx,Ny),
        levels,
        dx,
        dy)

    #check the sizes of inputs
    (size(mask) == (Nx,Ny)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to Nx x Ny (i.e. $Nx x $Ny)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent = c(Ny-1) ⊗ spI(Nx)
    ∂x = χ(Ny-2) ⊗ ∂1d(Nx-1,dx)
    ∂y = ∂1d(Ny-1,dy) ⊗ spI(Nx)
    dωt = wavelet_matrix(Ny,levels,"forward" ) ⊗ wavelet_matrix(Nx,levels,"forward")

    #fields stored on UGrid
    s = zeros(Nx,Ny)
    h = zeros(Nx,Ny)
    grounded_fraction = ones(Nx,Ny)
    βeff = zeros(Nx,Ny)
    dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)
    v = zeros(Nx,Ny)

    #size assertions
    @assert size(mask)==(Nx,Ny)
    n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == crop[mask[:],:]
    @assert spread == sparse(samp')
    @assert size(s)==(Nx,Ny)
    @assert size(h)==(Nx,Ny)
    @assert size(grounded_fraction)==(Nx,Ny)
    @assert size(βeff)==(Nx,Ny)
    @assert size(v)==(Nx,Ny)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return VGrid(
                Nx,
                Ny,
                mask,
                n,
                crop,
                samp,
                spread,
                cent,
                ∂x,
                ∂y,
                levels,
                dωt,
                s,
                h,
                grounded_fraction,
                βeff,
                dnegβeff,
                v)
end