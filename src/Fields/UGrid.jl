struct UGrid{T <: Real, N <: Integer}
                   nxu :: N                                    # Number of frid cells in x-direction in UGrid
                   nyu :: N                                    # Number of grid cells in y-direction in UGrid 
                 mask :: Array{Bool,2}                         # Mask specifying model domain wrt U grid 
           mask_inner :: Array{Bool,2}                         # Mask specifying interior of model domain wrt U grid 
            u_isfixed :: Array{Bool,2}                         # Mask specifying location of fixed u-velocity 
                    n :: N                                     # Total number of cells in model domain 
                   ni :: N                                     # Total number of cells in interior of model domain 
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
           samp_inner :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to interior of model domain 
               spread :: SparseMatrixCSC{T,N}                  # Spread matrix: take model domain to full domain
         spread_inner :: SparseMatrixCSC{T,N}                  # Spread matrix: take interior of model domain to full domain
                 cent :: KronType{T,N}                         # Map from U grid to H grid 
                centᵀ :: KronType{T,N}                         # Adjont of map from U grid to H grid 
                   ∂x :: KronType{T,N}                         # Matrix representation of differentiation wrt x 
                  ∂xᵀ :: KronType{T,N}                         # Adjoint of differentiation wrt x 
                   ∂y :: KronType{T,N}                         # Matrix representation of differentiation wrt y
                  ∂yᵀ :: KronType{T,N}                         # Adjoint of differentiation wrt y
               levels :: N                                     # Number of levels in the preconditioner
                    s :: Array{T,2}                            # Ice surface elevation
                    h :: Array{T,2}                            # Ice thickness
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction
                 βeff :: Array{T,2}                            # Effective β value on u grid(eqn 12 in Arthern 2015 JGeophysRes)
             dnegβeff :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (diagonal of βeff) 
                    u :: Array{T,2}                            # Ice velocities in x direction
             residual :: Array{T,2}                            # Residuals in x direction
                τsurf :: Array{T,2}                            # TauSurf Dirichlet calculated in inversion (see in Arthern 2015 JGeophysRes)
                   us :: Array{T,2}                               # Ice surface velocities in x direction
              end
    
"""
    UGrid(;
            nxu,
            nyu,
            mask = trues(nxu,nyu),
            u_isfixed = falses(nxu,nyu),
            u = zeros(nxu,nyu),
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
    - 'u_isfixed' Mask specifying where u velocities are fixed.
    - 'u' Values of u velocities (including fixed values).
    - levels: (required) Number of levels in the preconditioner
    - dx: (required) Grid spacing in the x direction
    - dy: (required) Grid spacing in the y direction
"""
function UGrid(;
                nxu,
                nyu,
                mask = trues(nxu,nyu),
                u_isfixed = falses(nxu,nyu),
                u = zeros(nxu,nyu),
                residual = zeros(nxu,nyu),
                levels,
                dx,
                dy,
                storage_only::Bool = false)

    #check the sizes of inputs
    (size(mask) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))
    (size(u_isfixed) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))
    (size(u) == (nxu,nyu)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxu x nyu (i.e. $nxu x $nyu)"))

    mask = convert(Array{Bool,2}, mask)
    u_isfixed = convert(Array{Bool,2}, u_isfixed)
    mask_inner = mask .& .! u_isfixed
    mask_inner = convert(Array{Bool,2}, mask_inner)

    #fields stored on UGrid
    s = zeros(nxu,nyu)
    h = zeros(nxu,nyu)
    grounded_fraction = ones(nxu,nyu)
    βeff = zeros(nxu,nyu)
    τsurf=zeros(nxu,nyu)
    us=zeros(nxu,nyu)

    if storage_only
        n = zero(nxu)
        ni = zero(nxu)
        crop = Diagonal(Float64[])
        samp = spzeros(Float64, 0, nxu * nyu)
        samp_inner = spzeros(Float64, 0, nxu * nyu)
        spread = spzeros(Float64, nxu * nyu, 0)
        spread_inner = spzeros(Float64, nxu * nyu, 0)
        kron = _storage_only_kron()
        cent = centᵀ = ∂x = ∂xᵀ = ∂y = ∂yᵀ = kron
        dnegβeff = Ref(Diagonal(Float64[]))
    else
        #construct operators
        n = count(mask)
        ni = count(mask_inner)
        crop = Diagonal(float(mask[:]))
        samp = sparse(1:n,(1:(nxu*nyu))[mask[:]],ones(n),n,nxu*nyu)
        samp_inner = sparse(1:ni,(1:(nxu*nyu))[mask_inner[:]],ones(ni),ni,nxu*nyu)
        spread = sparse(samp')
        spread_inner = sparse(samp_inner')
        cent =  spI(nyu) ⊗ c(nxu-1)
        centᵀ =  sparse(spI(nyu)') ⊗ sparse(c(nxu-1)')
        ∂x =  spI(nyu) ⊗ ∂1d(nxu-1,dx)
        ∂xᵀ =  sparse(spI(nyu)') ⊗ sparse(∂1d(nxu-1,dx)')
        ∂y =  ∂1d(nyu-1,dy) ⊗ χ(nxu-2)
        ∂yᵀ =  sparse(∂1d(nyu-1,dy)') ⊗ sparse(χ(nxu-2)')
        dnegβeff = Ref(crop*Diagonal(-βeff[:])*crop)

        #size assertions
        @assert n == count(mask)
        @assert ni == count(mask_inner)
    end

    @assert size(s)==(nxu,nyu)
    @assert size(h)==(nxu,nyu)
    @assert size(grounded_fraction)==(nxu,nyu)
    @assert size(βeff)==(nxu,nyu)
    @assert size(τsurf)==(nxu,nyu)
    @assert size(us)==(nxu,nyu)

    return UGrid(
                nxu,
                nyu,
                mask,
                mask_inner,
                u_isfixed,
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
                s,
                h,
                grounded_fraction,
                βeff,
                dnegβeff,
                u,
                residual,
                τsurf,
                us)
end