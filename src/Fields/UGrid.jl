struct UGrid{T <: Real, N <: Integer}
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
                    u :: Array{T,2}
end
    
"""
    UGrid(;
            Nx,
            Ny,
            mask = trues(Nx,Ny), 
            levels,
            dx,
            dy)

Construct a WAVI.jl UGrid with size (Nx,Ny)
UGrid stores fields that are defined on the problem's U grid. 
(Co-ordinates of HGrid stored in a Grid under xxhu yyu fields)

Keyword arguments
=================
    - 'Nx': (required) Number of grid cells in x-direction in UGrid (should be same as grid.nx + 1)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'Ny': (required) Number of grid cells in y-direction in HGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain with respect to U grid
    - levels: (required) Number of vertical levels 
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function UGrid(;
                Nx,
                Ny,
                mask = trues(Nx,Ny),
                levels,
                dx,
                dy)

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent =  spI(Ny) ⊗ c(Nx-1)
    ∂x = spI(Ny) ⊗ ∂1d(Nx-1,dx)
    ∂y =  ∂1d(Ny-1,dy) ⊗ χ(Nx-2)
    dωt = wavelet_matrix(Ny,levels,"forward" ) ⊗ wavelet_matrix(Nx,levels,"forward")

    #fields stored on UGrid
    s = zeros(Nx,Ny)
    h = zeros(Nx,Ny)
    grounded_fraction = ones(Nx,Ny)
    βeff = zeros(Nx,Ny)
    dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)
    u = zeros(Nx,Ny)


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
    @assert size(u)==(Nx,Ny)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return UGrid(
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
                u)
end