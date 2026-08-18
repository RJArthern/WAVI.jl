export UWavelets

#Struct to hold information on wavelet-grid (u-component).
struct UWavelets{T <: Real, N <: Integer} 
          nxuw :: N                                     # Number of grid points in x in UWavelets (equal to UGrid)
          nyuw :: N                                     # Number of grid points in y in UWavelets (equal to UGrid)
        mask :: Array{Bool,2}                           # Model domain on the U grid
           n :: Base.RefValue{N}                        # Number of grid points in domain
        crop :: Base.RefValue{Diagonal{T,Array{T,1}}}   # Crop matrix: diagonal matrix with mask entries on diag
        samp :: Base.RefValue{SparseMatrixCSC{T,N}}     # Sampling matrix: take full domain to model domain 
      spread :: Base.RefValue{SparseMatrixCSC{T,N}}     # Spread matrix: take model domain to full domain
      levels :: N                                       # Number of wavelet levels 
   index_map :: Matrix{Int}                             # Packed index of each kept wavelet coefficient; 0 if dropped
    wavelets :: Array{T,2}                              # Wavelet matrix
    correction_coarse::Base.RefValue{Vector{T} }        # cache to store coarse correction for multigrid                           
end

"""
    UWavelets(;
            nxuw,
            nyuw,
            mask = trues(nxuw,nyuw), 
            levels)

Construct a WAVI.jl UWavelets structure with size (nxuw,nyuw)
UWavelets stores wavelet matrix and associated fields for velocities in x direction

Keyword arguments
=================
    - 'nxuw': (required) Number of grid cells in x-direction in UWavelts 
    - 'nyuw': (required) Number of grid cells in y-direction in UWavelets 
    - 'mask': Mask specifying the model domain with respect to UWavelet grid (UGrid)
    - levels: (required) Number of levels in the preconditioner 
"""
function UWavelets(;   
                    nxuw,
                    nyuw,
                    mask = trues(nxuw,nyuw),
                    levels,
                    storage_only::Bool = false)
     
    (size(mask) == (nxuw,nyuw)) || throw(DimensionMismatch("Sizes of inputs to UWavelets must all be equal to nxuw x nyuw (i.e. $nxuw x $nyuw)"))

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    wavelets = zeros(nxuw,nyuw); @assert size(wavelets)==(nxuw,nyuw)
    index_map = zeros(Int, nxuw, nyuw)

    if storage_only
        n = Ref(0)
        crop = Ref(Diagonal(Float64[]))
        samp = Ref(spzeros(Float64, 0, nxuw * nyuw))
        spread = Ref(spzeros(Float64, nxuw * nyuw, 0))
        correction_coarse = Ref(Float64[])
    else
        #compute non-inputs
        n = Ref(count(mask));  @assert n[] == count(mask)
        crop = Ref(Diagonal(float(mask[:])))
        samp  = Ref(sparse(1:n[],(1:(nxuw*nyuw))[mask[:]],ones(n[]),n[],nxuw*nyuw))
        spread = Ref(sparse(samp[]'))
        fill_index_map!(index_map, mask)
        correction_coarse = Ref(zeros(n[])); @assert length(correction_coarse[]) == n[]
    end

    return UWavelets(
                    nxuw,
                    nyuw,
                    mask,
                    n,
                    crop, 
                    samp, 
                    spread, 
                    levels,
                    index_map,
                    wavelets,
                    correction_coarse)
end
