struct CGrid{T <: Real, N <: Integer}
        nxc :: N
        nyc :: N
      mask :: Array{Bool,2} 
         n :: N 
      crop :: Diagonal{T,Array{T,1}}
      samp :: SparseMatrixCSC{T,N} 
    spread :: SparseMatrixCSC{T,N}
      cent :: KronType{T,N}
  dneghηav :: Base.RefValue{Diagonal{T,Array{T,1}}}    # Rheological operator (-h × ηav)
end

function CGrid(;
                nxc,
                nyc,
                mask = trues(nxc,nyc))

    #check the sizes of inputs
    (size(mask) == (nxc,nyc)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to nxc x nyc (i.e. $nxc x $nyc)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent = sparse(c(nyc)') ⊗ sparse(c(nxc)')
    dneghηav = Ref(crop*Diagonal(zeros(nxc*nyc))*crop)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

        
"""
CGrid(;
        nxc,
        nyc,
        mask = trues(nxc,nyc))

Construct a WAVI.jl CGrid with size (nxc,nyc)
CGrid stores operators required for the CGrid. 
(Co-ordinates of CGrid stored in a Grid under xxc and yyc fields)

Keyword arguments
=================
- 'nxc': (required) Number of grid cells in x-direction in CGrid (should be same as grid.nx -1)
        Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
- 'nyc': (required) Number of grid cells in y-direction in CGrid (should be same as grid.ny - 1)
- 'mask': Mask specifying the model domain with respect to C grid
"""
return CGrid(
            nxc,
            nyc, 
            mask,
            n,
            crop, 
            samp, 
            spread,
            cent,
            dneghηav)
end
    
