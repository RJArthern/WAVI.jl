
#Type to allow non-allocating multiplication by Kronecker products

@with_kw struct KroneckerProduct{T <: Real,TA <: AbstractMatrix{T},TB <: AbstractMatrix{T}} <: LinearMap{T}
    A::TA
    B::TB
    temp :: Matrix{T} 
end

function KroneckerProduct(A,B) 
    if size(B,2)*size(A,1) <= size(B,1)*size(A,2)
        return KroneckerProduct(A,B,zeros(size(B,2),size(A,1)))
    else
        return KroneckerProduct(A,B,zeros(size(B,1),size(A,2)))
    end
end

Base.size(K::KroneckerProduct) = size(K.A,1)*size(K.B,1), size(K.A,2)*size(K.B,2)
Base.size(K::KroneckerProduct,i::Int) = size(K.A,i)*size(K.B,i)

function Base.:(*)(K::KroneckerProduct,x::AbstractVector)
    @assert size(K,2)==size(x,1)
    y=zeros(promote_type(eltype(K),eltype(x)),size(K,1))
    mul!(y,K,x)
end

function mul!(y::AbstractVecOrMat,K::KroneckerProduct,x::AbstractVector) 
    @unpack A, B, temp = K
    (size(y,1)==size(K,1) && size(x,1)==size(K,2) && size(y,2) == size(x,2))||
         throw(DimensionMismatch("Inconsistent dimensions of Kronecker product."))
    ma, na = size(A)
    mb, nb = size(B)
    X = Base.ReshapedArray(x, (nb, na),())
    Y = Base.ReshapedArray(y, (mb, ma),())
    if nb*ma <= mb*na
        mul!(temp, X, transpose(A))
        mul!(Y, B, temp)
    else
        mul!(temp, B, X)
        mul!(Y, temp, transpose(A))
    end
    return y
end

⊗(A::SparseMatrixCSC,B::SparseMatrixCSC) = KroneckerProduct(A,sparse(B')') 
const KronType{T,N} = KroneckerProduct{T, SparseMatrixCSC{T, N}, Adjoint{T,SparseMatrixCSC{T, N}}} where {T,N}

#⊗(A::SparseMatrixCSC,B::SparseMatrixCSC) = KroneckerProduct(ThreadedSparseMatrixCSC(A),ThreadedSparseMatrixCSC(B)) 
#Type alias just for abreviation
#const KronType{T,N} = KroneckerProduct{T, ThreadedSparseMatrixCSC{T, N, SparseMatrixCSC{T, N}}, 
#  ThreadedSparseMatrixCSC{T, N, SparseMatrixCSC{T, N}}} where {T,N}