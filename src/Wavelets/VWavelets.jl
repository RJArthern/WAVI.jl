export VWavelets

#Struct to hold information on wavelet-grid (v-component).
struct VWavelets{T <: Real, N <: Integer}
        nxvw :: N                                       # Number of grid points in x in UWavelets (equal to UGrid)
        nyvw :: N                                       # Number of grid points in y in UWavelets (equal to UGrid)
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
    VWavelets(;
            nxvw,
            nyvw,
            mask = trues(nxvw,nyvw), 
            levels)

Construct a WAVI.jl VWavelets structure with size (nxvw,nyvw)
UWavelets stores wavelet matrix and associated fields for velocities in y direction

Keyword arguments
=================
    - 'nxvw': (required) Number of grid cells in x-direction in VWavelets 
    - 'nyvw': (required) Number of grid cells in y-direction in VWavelets 
    - 'mask': Mask specifying the model domain with respect to UWavelet grid (VGrid)
    - levels: (required) Number of levels in the preconditioner 
"""
function VWavelets(;   
                    nxvw,
                    nyvw,
                    mask = trues(nxvw,nyvw),
                    levels,
                    storage_only::Bool = false)
     
    (size(mask) == (nxvw,nyvw)) || throw(DimensionMismatch("Sizes of inputs to UWavelets must all be equal to nxvw x nyvw (i.e. $nxvw x $nyvw)"))

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    wavelets = zeros(nxvw,nyvw); @assert size(wavelets)==(nxvw,nyvw)
    index_map = zeros(Int, nxvw, nyvw)

    if storage_only
        n = Ref(0)
        crop = Ref(Diagonal(Float64[]))
        samp = Ref(spzeros(Float64, 0, nxvw * nyvw))
        spread = Ref(spzeros(Float64, nxvw * nyvw, 0))
        correction_coarse = Ref(Float64[])
    else
        #compute non-inputs
        n = Ref(count(mask));  @assert n[] == count(mask)
        crop = Ref(Diagonal(float(mask[:])))
        samp  = Ref(sparse(1:n[],(1:(nxvw*nyvw))[mask[:]],ones(n[]),n[],nxvw*nyvw))
        spread = Ref(sparse(samp[]'))
        fill_index_map!(index_map, mask)
        correction_coarse = Ref(zeros(n[])); @assert length(correction_coarse[])==n[]
    end

    return VWavelets(
                    nxvw,
                    nyvw,
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
