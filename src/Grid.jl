struct Grid{T <: Real, N <: Integer} <: AbstractGrid{T,N}
                    nx :: N             # Number of x gridpoints
                    ny :: N             # Number of y gridpoints
                    nσ :: N             # Number of levels in the vertical
                    dx :: T             # Grid spacing in x
                    dy :: T             # Grid spacing in y
                    x0 :: T             # X co-ordinate of grid origin
                    y0 :: T             # Y co-ordinate of grid origin
                h_mask :: Array{Bool,2} # Mask defining domain points within grid
              u_iszero :: Array{Bool,2} # Locations of zero u velocity points 
              v_iszero :: Array{Bool,2} # Locations of zero v velocity points
                   xxh :: Array{T,2}    # x co-ordinates matrix of h grid
                   yyh :: Array{T,2}    # y co-ordinates matrix of h grid
                   xxu :: Array{T,2}    # x co-ordinates matrix of u grid
                   yyu :: Array{T,2}    # y co-ordinates matrix of u grid
                   xxv :: Array{T,2}    # x co-ordinates matrix of v grid
                   yyv :: Array{T,2}    # y co-ordinates matrix of v grid
                   xxc :: Array{T,2}    # x co-ordinates matrix of c grid
                   yyc :: Array{T,2}    # y co-ordinates matrix of c grid
                     σ :: Vector{T}     # Dimensionless levels in the vertical
                     ζ :: Vector{T}     # Reverse dimensionless levels in the vertical
    quadrature_weights :: Vector{T}     # Quadrature weights for integration
                   Cxl :: N              #lower x extent value of coupled child domain
                   Cxu :: N              #upper x extent value of coupled child domain
                   Cyl :: N              #lower y extent value of coupled child domain
                   Cyu :: N              #upper y extent value of coupled child domain
end

"""
    Grid(; 
    nx = 80,
    ny = 10,
    dx = 8000.0,
    dy = 8000.0,
    nσ = 4,
    x0 = 0.0,
    y0 = -40000.0,
    h_mask = trues(nx,ny),
    u_iszero = falses(nx+1,ny),
    v_iszero = falses(nx,ny+1),
    Cxl = 1,
    Cxu = nx,
    Cyl = 1,
    Cyu = ny)

Construct a WAVI.jl grid.

Keyword arguments
=================

    - 'nx': number of x grid points
    - 'ny': number of y grid points
    - 'dx': grid spacing in x 
    - 'dy': grid spacing in y 
    - 'nσ': number of levels in the vertical
    - 'x0': grid origin x co-ordinate 
    - 'y0': grid origin y co-ordinate
    - 'h_mask': Mask defining domain points within grid
    - 'u_iszero': Locations of zero u velocity points
    - 'v_iszero': Locations of zero v velocity points
"""

#grid constructor 
function Grid(; 
    nx = 80,
    ny = 10,
    dx = 8000.0,
    dy = 8000.0,
    nσ = 4,
    x0 = 0.0,
    y0 = -40000.0,
    h_mask = trues(nx,ny),
    u_iszero = falses(nx+1,ny),
    v_iszero = falses(nx,ny+1),
    Cxl = 1,
    Cxu = 80,
    Cyl = 1,
    Cyu = 10)

#check the sizes of inputs
@assert size(h_mask)==(nx,ny);@assert h_mask == clip(h_mask)
@assert size(u_iszero)==(nx+1,ny)
@assert size(v_iszero)==(nx,ny+1)

#map bit arrays to boolean
h_mask = convert(Array{Bool,2}, h_mask)
u_iszero = convert(Array{Bool,2}, u_iszero)
v_iszero = convert(Array{Bool,2}, v_iszero)


#compute grid co-ordinates
xxh=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xxh)==(nx,ny)
yyh=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yyh)==(nx,ny)

xxu=[x0+(i-1.0)*dx for i=1:nx, j=1:ny]; @assert size(xxu)==(nx,ny)
yyu=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yyu)==(nx,ny)

xxv=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xxv)==(nx,ny)
yyv=[y0+(j-1.0)*dy for i=1:nx, j=1:ny]; @assert size(yyv)==(nx,ny)

xxc=[x0+i*dx for i=1:nx, j=1:ny]; @assert size(xxc)==(nx,ny)
yyc=[y0+j*dy for i=1:nx, j=1:ny]; @assert size(yyc)==(nx,ny)

#sigma grid info
σ = collect(range(0.0,length=nσ,stop=1.0)); @assert length(σ) == nσ
ζ = one(eltype(σ)) .- σ ; @assert length(ζ) == nσ
quadrature_weights = [0.5;ones(nσ-2);0.5]/(nσ-1); @assert length(quadrature_weights) == nσ

return Grid(nx,
            ny,
            nσ,
            dx,
            dy,
            x0,
            y0,
            h_mask,
            u_iszero,
            v_iszero,
            xxh,
            yyh,
            xxu,
            yyu,
            xxv,
            yyv,
            xxc,
            yyc,
            σ,
            ζ,
            quadrature_weights)
end
