struct UGrid{T <: Real, N <: Integer}
    Nx::N
    Ny::N
    mask::Array{Bool,2}
    n::N 
    crop::Diagonal{T,Array{T,1}}
    samp::SparseMatrixCSC{T,N}
    spread::SparseMatrixCSC{T,N}
    cent::KronType{T,N}
    ∂x::KronType{T,N}
    ∂y::KronType{T,N}
    levels::N
    dwt::KronType{T,N}
    s::Array{T,2} 
    h::Array{T,2}
    grounded_fraction::Array{T,2}
    βeff::Array{T,2}
    dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} 
    u::Array{T,2}
end
    
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