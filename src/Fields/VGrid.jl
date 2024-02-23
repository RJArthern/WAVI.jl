struct VGrid{T <: Real, N <: Int}
                  nxv :: N                                     # Number of frid cells in x-direction in UGrid
                  nyv :: N                                     # Number of grid cells in y-direction in UGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt V grid 
           mask_inner :: Array{Bool,2}                         # Mask specifying interior of model domain wrt V grid 
            v_isfixed :: Array{Bool,2}                         # Mask specifying location of fixed v-velocity 
                    n :: N                                     # Total number of cells in model domain 
                   ni :: N                                     # Total number of cells in interior of model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
           samp_inner :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to interior of model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
         spread_inner :: SparseMatrixCSC{T,N}                  # Spread matrix: take interior of model domain to full domain
                 cent :: KronType{T,N}                         # Map from V grid to H grid 
                centᵀ :: KronType{T,N}                         # Adjoint of map from V grid to H grid 
                   ∂x :: KronType{T,N}                         # Matrix representation of differentiation wrt x 
                  ∂xᵀ :: KronType{T,N}                         # Adjoint of differentiation wrt x 
                   ∂y :: KronType{T,N}                         # Matrix representation of differentiation wrt y
                  ∂yᵀ :: KronType{T,N}                         # Adjoint of differentiation wrt y
               levels :: N                                     # Number of levels in the preconditioner
                  dwt :: KronType{T,N}                         # Wavelet matrix product on u grid 
                    s :: Array{T,2}                            # Ice surface elevation
                    h :: Array{T,2}                            # Ice thickness
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction
                 βeff :: Array{T,2}                            # Effective β value on u grid(eqn 12 in Arthern 2015 JGeophysRes)
             dnegβeff :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (diagonal of βeff) 
                    v :: Array{T,2}                            # Ice velocity in the y-direction
             residual :: Array{T,2}                            # Residuals in the y-direction
            end

    
"""
    VGrid(;
            nxv,
            nyv,
            mask = trues(nxv,nyv), 
            v_isfixed = falses(nxv,nyv),
            v = zeros(nxv,nyv),
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
    - 'v_isfixed' Mask specifying where v velocities are fixed
    - 'v' Values of v velocities (including fixed values).
    - levels: (required) Number of levels in the preconditioner 
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function VGrid(;
        nxv,
        nyv,
        mask = trues(nxv,nyv),
        v_isfixed = falses(nxv,nyv),
        v = zeros(nxv,nyv),
        residual = zeros(nxv,nyv),
        levels,
        dx,
        dy)

    #check the sizes of inputs
    (size(mask) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to VGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))
    (size(v_isfixed) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to VGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))
    (size(v) == (nxv,nyv)) || throw(DimensionMismatch("Sizes of inputs to VGrid must all be equal to nxv x nyv (i.e. $nxv x $nyv)"))

    #construct operators
    n = count(mask)
    mask_inner = mask .& .! v_isfixed
    ni = count(mask_inner)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxv*nyv))[mask[:]],ones(n),n,nxv*nyv)
    samp_inner = sparse(1:ni,(1:(nxv*nyv))[mask_inner[:]],ones(ni),ni,nxv*nyv)
    spread = sparse(samp')
    spread_inner = sparse(samp_inner')
    cent = c(nyv-1) ⊗ spI(nxv)
    centᵀ = sparse(c(nyv-1)') ⊗ sparse(spI(nxv)')
    ∂x = χ(nyv-2) ⊗ ∂1d(nxv-1,dx)
    ∂xᵀ = sparse(χ(nyv-2)') ⊗ sparse(∂1d(nxv-1,dx)')
    ∂y = ∂1d(nyv-1,dy) ⊗ spI(nxv)
    ∂yᵀ = sparse(∂1d(nyv-1,dy)') ⊗ sparse(spI(nxv)')
    dωt = wavelet_matrix(nyv,levels,"forward" ) ⊗ wavelet_matrix(nxv,levels,"forward")

    #fields stored on UGrid
    s = zeros(nxv,nyv)
    h = zeros(nxv,nyv)
    grounded_fraction = ones(nxv,nyv)
    βeff = zeros(nxv,nyv)
    dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)

    #size assertions
    @assert n == count(mask)
    @assert ni == count(mask_inner)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxv*nyv))[mask[:]],ones(n),n,nxv*nyv)
    @assert samp_inner == sparse(1:ni,(1:(nxv*nyv))[mask_inner[:]],ones(ni),ni,nxv*nyv)
    @assert spread == sparse(samp')
    @assert spread_inner == sparse(samp_inner')
    @assert size(s)==(nxv,nyv)
    @assert size(h)==(nxv,nyv)
    @assert size(grounded_fraction)==(nxv,nyv)
    @assert size(βeff)==(nxv,nyv)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    mask_inner = convert(Array{Bool,2}, mask_inner)
    v_isfixed = convert(Array{Bool,2}, v_isfixed)


    return VGrid(
                nxv,
                nyv,
                mask,
                mask_inner,
                v_isfixed,
                n,
                ni,
                crop,
                samp,
                samp_inner,
                spread,
                spread_inner,
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
                v,
                residual)
end
