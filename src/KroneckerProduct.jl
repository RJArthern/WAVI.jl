
#Type to allow non-allocating multiplication by Kronecker products

"""
  KroneckerProduct

Struct to represent a lazy representation of the Kronecker product K = A ⊗ B. 
The matrix is never constructed, but K can multiply a vector.
The struct contains a temporary cache used in multiplication. 

"""
@with_kw struct KroneckerProduct{T <: Real,TA <: AbstractMatrix{T},TB <: AbstractMatrix{T}} <: LinearMap{T}
    A::TA
    B::TB
    temp :: Matrix{T} 
end


"""
  KroneckerProduct(A,B)

Constructor for a lazy representation of the Kronecker product K = A ⊗ B. 
The matrix is never constructed, but K can multiply a vector. 

"""
function KroneckerProduct(A,B) 
    if size(B,2)*size(A,1) <= size(B,1)*size(A,2)
        return KroneckerProduct(A,B,zeros(size(B,2),size(A,1)))
    else
        return KroneckerProduct(A,B,zeros(size(B,1),size(A,2)))
    end
end

Base.size(K::KroneckerProduct) = size(K.A,1)*size(K.B,1), size(K.A,2)*size(K.B,2)
Base.size(K::KroneckerProduct,i::Int) = size(K.A,i)*size(K.B,i)


"""
    Base.:(*)(K::KroneckerProduct,x::AbstractVector)

Function to multiply a vector by a Kronecker product

"""
function Base.:(*)(K::KroneckerProduct,x::AbstractVector)
    @assert size(K,2)==size(x,1)
    y=zeros(promote_type(eltype(K),eltype(x)),size(K,1))
    mul!(y,K,x)
end


"""
    mul!(y::AbstractVecOrMat,K::KroneckerProduct,x::AbstractVector) 

In-place function to multiply a vector by a Kronecker product

"""
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

"""
⊗(A::SparseMatrixCSC,B::SparseMatrixCSC)

Constructs a lazy representation of the Kronecker product K = A ⊗ B. 
The matrix is never constructed, but K can multiply a vector. 

"""
⊗(A::SparseMatrixCSC,B::SparseMatrixCSC) = KroneckerProduct(A,sparse(B')') 

#Type alias for abreviation
const KronType{T,N} = KroneckerProduct{T, SparseMatrixCSC{T, N}, Adjoint{T,SparseMatrixCSC{T, N}}} where {T,N}

