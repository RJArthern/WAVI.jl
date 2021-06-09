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
    h_mask = nothing,
    u_iszero = nothing,
    v_iszero = nothing)

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
    - 'quadrature_weights': weights associated with sigma levels used in quadrature scheme
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
    h_mask = nothing,
    u_iszero = nothing,
    v_iszero = nothing,
    quadrature_weights = nothing)

#check integer inputs
((typeof(nx) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in x direction (nx) must a positive integer larger than one")) 
((typeof(ny) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in y direction (ny)  must a positive integer larger than one")) 
((typeof(nσ) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in vertical (nσ)  must a positive integer larger than one")) 

#if nx, ny, nσ  of correct type, assemble h_mask, u_iszero, v_iszero (if not passed)
(~(h_mask === nothing)) || (h_mask = trues(nx,ny))
(~(u_iszero === nothing))|| (u_iszero = falses(nx+1,ny))
(~(v_iszero === nothing)) || (v_iszero = falses(nx, ny+1))
(~(quadrature_weights === nothing) || (quadrature_weights = [0.5;ones(nσ-2);0.5]/(nσ-1)))
 
#check the sizes of inputs
size(h_mask)==(nx,ny) || throw(DimensionMismatch("h_mask size must be (nx x ny) (i.e. $nx x $ny)"))
size(quadrature_weights) == (nσ,) || throw(DimensionMismatch("Input quadrate weighs are size $size(quadrature_weights). quadrature weights must have size (nσ,) (i.e. ($nσ,))"))
size(u_iszero)==(nx+1,ny) || throw(DimensionMismatch("u_iszero size must be size of U grid (nx+1 x ny) (i.e. $(nx+1) x $ny)"))
size(v_iszero)==(nx,ny+1) || throw(DimensionMismatch("v_iszero size must be size of V grid (nx x ny+1) (i.e. $nx x $(ny+1)"))

#map bit arrays to boolean
try
    h_mask = convert(Array{Bool,2}, h_mask)
    @assert h_mask == clip(h_mask)
catch 
    throw(ArgumentError("h_mask must be Boolean (or equivalent)"))
end

try 
    u_iszero = convert(Array{Bool,2}, u_iszero)
catch 
    throw(ArgumentError("u_iszero must be Boolean (or equivalent)"))
end

try
    v_iszero = convert(Array{Bool,2}, v_iszero)
catch
    throw(ArgumentError("v_iszero must be Boolean (or equivalent)"))
end

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
