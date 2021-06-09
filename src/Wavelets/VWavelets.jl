#Struct to hold information on wavelet-grid (v-component).
struct VWavelets{T <: Real, N <: Integer}
          Nx :: N                                       # Number of grid points in x in UWavelets (equal to UGrid)
          Ny :: N                                       # Number of grid points in y in UWavelets (equal to UGrid)
        mask :: Array{Bool,2}                           # Model domain on the U grid
           n :: Base.RefValue{N}                        # Number of grid points in domain
        crop :: Base.RefValue{Diagonal{T,Array{T,1}}}   # Crop matrix: diagonal matrix with mask entries on diag
        samp :: Base.RefValue{SparseMatrixCSC{T,N}}     # Sampling matrix: take full domain to model domain 
      spread :: Base.RefValue{SparseMatrixCSC{T,N}}     # Spread matrix: take model domain to full domain
      levels :: N                                       # Number of wavelet levels 
        idwt :: KronType{T,N}                           # Wavelet matrix cross produce
    wavelets :: Array{T,2}                              # Wavelet matrix
end


"""
    VWavelets(;
            Nx,
            Ny,
            mask = trues(Nx,Ny), 
            levels)

Construct a WAVI.jl VWavelets structure with size (Nx,Ny)
UWavelets stores wavelet matrix and associated fields for velocities in y direction

Keyword arguments
=================
    - 'Nx': (required) Number of grid cells in x-direction in VWavelets 
    - 'Ny': (required) Number of grid cells in y-direction in VWavelets 
    - 'mask': Mask specifying the model domain with respect to UWavelet grid (VGrid)
    - levels: (required) Number of levels in the preconditioner 
"""
function VWavelets(;   
                    Nx,
                    Ny,
                    mask = trues(Nx,Ny),
                    levels)
     
    (size(mask) == (Nx,Ny)) || throw(DimensionMismatch("Sizes of inputs to UWavelets must all be equal to Nx x Ny (i.e. $Nx x $Ny)"))

    #compute non-inputs 
    n = Ref(count(mask));  @assert n[] == count(mask)
    crop = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]));
    samp  = Ref(crop[][mask[:],:]); @assert samp[] == crop[][mask[:],:]
    spread = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
    idwt =  wavelet_matrix(Ny,levels,"reverse" ) ⊗ wavelet_matrix(Nx,levels,"reverse")
    wavelets = zeros(Nx,Ny); @assert size(wavelets)==(Nx,Ny)
    
    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)

    return VWavelets(
                    Nx,
                    Ny,
                    mask,
                    n,
                    crop, 
                    samp, 
                    spread, 
                    levels,
                    idwt,
                    wavelets)
end
    