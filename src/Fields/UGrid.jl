struct UGrid{T <: Real, N <: Integer, K1 <: KronType{T},
                                      K2 <: KronType{T},
                                      K3 <: KronType{T},
                                      K4 <: KronType{T},
                                      K5 <: KronType{T},
                                      K6 <: KronType{T},
                                      K7 <: KronType{T}}
                   nxu :: N                                    # Number of frid cells in x-direction in UGrid
                   nyu :: N                                    # Number of grid cells in y-direction in UGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
                    n :: N                                     # Total number of cells in model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
                 cent :: K1                                    # Map from U grid to H grid 
                centᵀ :: K2                                    # Adjont of map from U grid to H grid 
                   ∂x :: K3                                    # Matrix representation of differentiation wrt x 
                  ∂xᵀ :: K4                                    # Adjoint of differentiation wrt x 
                   ∂y :: K5                                    # Matrix representation of differentiation wrt y
                  ∂yᵀ :: K6                                    # Adjoint of differentiation wrt y
               levels :: N                                     # Number of levels in the preconditioner
                  dwt :: K7                                    # Wavelet matrix product on u grid 
                    s :: Array{T,2}                            # Ice surface elevation
                    h :: Array{T,2}                            # Ice thickness
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction
                 βeff :: Array{T,2}                            # Effective β value on u grid(eqn 12 in Arthern 2015 JGeophysRes)
             dnegβeff :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (diagonal of βeff) 
                    u :: Array{T,2}                            # Ice velocities in x direction
end
    
"""
    UGrid(;
            nxu,
            nyu,
            mask = trues(nxu,nyu), 
            levels,
            dx,
            dy)

Construct a WAVI.jl UGrid with size (nxu,nyu)
UGrid stores fields that are defined on the problem's U grid. 
(Co-ordinates of UGrid stored in a Grid under xxu and yyu fields)

Keyword arguments
=================
    - 'nxu': (required) Number of grid cells in x-direction in UGrid (should be same as grid.nx + 1)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyu': (required) Number of grid cells in y-direction in UGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain with respect to U grid
    - levels: (required) Number of levels in the preconditioner
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function UGrid(;
                nxu,
                nyu,
                mask = trues(nxu,nyu),
                levels,
                dx,
                dy)

    #check the sizes of inputs
    (size(mask) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent =  spI(nyu) ⊗ c(nxu-1)
    centᵀ =  sparse(spI(nyu)') ⊗ sparse(c(nxu-1)')
    ∂x =  spI(nyu) ⊗ ∂1d(nxu-1,dx)
    ∂xᵀ =  sparse(spI(nyu)') ⊗ sparse(∂1d(nxu-1,dx)')
    ∂y =  ∂1d(nyu-1,dy) ⊗ χ(nxu-2)
    ∂yᵀ =  sparse(∂1d(nyu-1,dy)') ⊗ sparse(χ(nxu-2)')
    dωt = wavelet_matrix(nyu,levels,"forward" ) ⊗ wavelet_matrix(nxu,levels,"forward")

    #fields stored on UGrid
    s = zeros(nxu,nyu)
    h = zeros(nxu,nyu)
    grounded_fraction = ones(nxu,nyu)
    βeff = zeros(nxu,nyu)
    dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)
    u = zeros(nxu,nyu)


    #size assertions
    @assert size(mask)==(nxu,nyu)
    n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == crop[mask[:],:]
    @assert spread == sparse(samp')
    @assert size(s)==(nxu,nyu)
    @assert size(h)==(nxu,nyu)
    @assert size(grounded_fraction)==(nxu,nyu)
    @assert size(βeff)==(nxu,nyu)
    @assert size(u)==(nxu,nyu)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return UGrid(
                nxu,
                nyu,
                mask,
                n,
                crop,
                samp,
                spread,
                cent,
                centᵀ,
                ∂x,
                ∂xᵀ,
                ∂y,
                ∂yᵀ,
                levels,
                dωt,
                s,
                h,
                grounded_fraction,
                βeff,
                dnegβeff,
                u)
end