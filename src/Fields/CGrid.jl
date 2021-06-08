struct CGrid{T <: Real, N <: Integer}
        Nx :: N
        Ny :: N
      mask :: Array{Bool,2} 
         n :: N 
      crop :: Diagonal{T,Array{T,1}}
      samp :: SparseMatrixCSC{T,N} 
    spread :: SparseMatrixCSC{T,N}
      cent :: KronType{T,N}
end

function CGrid(;
                Nx,
                Ny,
                mask = trues(Nx,Ny))

    #check the sizes of inputs
    (size(mask) == (Nx,Ny)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to Nx x Ny (i.e. $Nx x $Ny)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = crop[mask[:],:]
    spread = sparse(samp')
    cent = sparse(c(Ny)') ⊗ sparse(c(Nx)')

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

        
"""
CGrid(;
        Nx,
        Ny,
        mask = trues(Nx,Ny))

Construct a WAVI.jl CGrid with size (Nx,Ny)
CGrid stores operators required for the CGrid. 
(Co-ordinates of CGrid stored in a Grid under xxc and yyc fields)

Keyword arguments
=================
- 'Nx': (required) Number of grid cells in x-direction in CGrid (should be same as grid.nx -1)
        Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
- 'Ny': (required) Number of grid cells in y-direction in CGrid (should be same as grid.ny - 1)
- 'mask': Mask specifying the model domain with respect to C grid
"""
return CGrid(
            Nx,
            Ny, 
            mask,
            n,
            crop, 
            samp, 
            spread,
            cent)
end
    
