struct VGrid{T <: Real, N <: Int}
                   nxv :: N                                     # Number of frid cells in x-direction in UGrid
                   nyv :: N                                     # Number of grid cells in y-direction in UGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
                    n :: N                                     # Total number of cells in model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
                 cent :: KronType{T,N}                         # Map from U grid to H grid 
                   ∂x :: KronType{T,N}                         # Matrix representation of differentiation wrt x 
                   ∂y :: KronType{T,N}                         # Matrix representation of differentiation wrt y
               levels :: N                                     # Number of levels in the preconditioner
                  dwt :: KronType{T,N}                         # Wavelet matrix product on u grid 
                    s :: Array{T,2}                            # Ice surface elevation
                    h :: Array{T,2}                            # Ice thickness
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction
                 βeff :: Array{T,2}                            # Effective β value on u grid(eqn 12 in Arthern 2015 JGeophysRes)
             dnegβeff :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (diagonal of βeff) 
                    v :: Array{T,2}                            # Ice velocity in the y-direction
end

    
"""
    VGrid(;
            nxv,
            nyv,
            mask = trues(nxv,nyv), 
            levels,
            dx,
            dy)

Construct a WAVI.jl VGrid with size (nxv,nyv)
VGrid stores fields that are defined on the problem's V grid. 
(Co-ordinates of HGrid stored in a Grid under xxv and yyv fields)

Keyword arguments
=================
    - 'nxv': (required) Number of grid cells in x-direction in VGrid (should be same as grid.nx)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyv': (required) Number of grid cells in y-direction in VGrid (should be same as grid.ny + 1)
    - 'mask': Mask specifying the model domain with respect to V grid
    - levels: (required) Number of levels in the preconditioner 
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function VGrid(;
        nxv,
        nyv,
        mask = trues(nxv,nyv),
        levels,
        dx,
        dy)

    #check the sizes of inputs
    (size(mask) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent = c(nyv-1) ⊗ spI(nxv)
    ∂x = χ(nyv-2) ⊗ ∂1d(nxv-1,dx)
    ∂y = ∂1d(nyv-1,dy) ⊗ spI(nxv)
    dωt = wavelet_matrix(nyv,levels,"forward" ) ⊗ wavelet_matrix(nxv,levels,"forward")

    #fields stored on UGrid
    s = zeros(nxv,nyv)
    h = zeros(nxv,nyv)
    grounded_fraction = ones(nxv,nyv)
    βeff = zeros(nxv,nyv)
    dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)
    v = zeros(nxv,nyv)

    #size assertions
    @assert size(mask)==(nxv,nyv)
    n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == crop[mask[:],:]
    @assert spread == sparse(samp')
    @assert size(s)==(nxv,nyv)
    @assert size(h)==(nxv,nyv)
    @assert size(grounded_fraction)==(nxv,nyv)
    @assert size(βeff)==(nxv,nyv)
    @assert size(v)==(nxv,nyv)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return VGrid(
                nxv,
                nyv,
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