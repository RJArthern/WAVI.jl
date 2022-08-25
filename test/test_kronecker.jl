using WAVI, Test, SparseArrays, LinearAlgebra

@testset "KroneckerProduct" begin
    @info "Testing Kronecker product ...."

    nbig=1000
    nsmall=100
    sparsity=0.01
    A=sprandn(nbig,nsmall,sparsity)
    B=sprandn(nsmall,nbig,sparsity)
    x=randn(nbig*nsmall)

    K0=kron(A,B)
    K1=WAVI.KroneckerProduct(A,B)
    y0=zeros(nbig*nsmall)
    y1=zeros(nbig*nsmall)
    mul!(y0,K0,x)
    mul!(y1,K1,x)
    @test y0 ≈ y1 

    K0=kron(B,A)
    K1=WAVI.KroneckerProduct(B,A)
    y0=zeros(nbig*nsmall)
    y1=zeros(nbig*nsmall)
    mul!(y0,K0,x)
    mul!(y1,K1,x)
    @test y0 ≈ y1 

    xwrong=randn(nbig*nsmall+1)
    @test_throws DimensionMismatch mul!(y1,K1,xwrong)
    ywrong=zeros(nbig*nsmall+1)
    @test_throws DimensionMismatch mul!(ywrong,K1,x)
    ywrong=zeros(nbig*nsmall,2)
    @test_throws DimensionMismatch mul!(ywrong,K1,x)

end