module Grids

import Base: show, size
export Grid,
    reconstruct_on_grid,
    reconstruct_on_subdomain,
    spatial_on_subdomain,
    field_on_grid,
    subdomain_index_ranges,
    grid_index_ranges,
    forcing_index_ranges

using WAVI: AbstractGrid

struct Grid{T <: Real, N <: Integer} <: AbstractGrid{T,N}
                    nx :: N             # Number of x gridpoints
                    ny :: N             # Number of y gridpoints
                    nσ :: N             # Number of levels in the vertical
                    dx :: T             # Grid spacing in x
                    dy :: T             # Grid spacing in y
                    x0 :: T             # X co-ordinate of grid origin
                    y0 :: T             # Y co-ordinate of grid origin
                h_mask :: Array{Bool,2} # Mask defining domain points within grid
             h_isfixed :: Array{Bool,2} # Mask defining locations of fixed thickness within grid
 hyd_potential_isfixed :: Array{Bool,2} # Mask defining locations of fixed hydraulic potential at the bed within grid
              u_iszero :: Array{Bool,2} # Locations of zero u velocity points 
              v_iszero :: Array{Bool,2} # Locations of zero v velocity points
             u_isfixed :: Array{Bool,2} # Locations of fixed u velocity points 
             v_isfixed :: Array{Bool,2} # Locations of fixed v velocity points
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
              basin_ID :: Array{T,2}    # grid of IDs for different basins
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
    h_isfixed = nothing,
    hyd_potential_isfixed = nothing,
    u_iszero = nothing,
    v_iszero = nothing,
    basin_ID = nothing,
    u_isfixed = nothing,
    v_isfixed = nothing)

Construct a WAVI.jl grid.

Keyword arguments
=================
- `nx`: number of x grid points
- `ny`: number of y grid points
- `dx`: grid spacing in x 
- `dy`: grid spacing in y 
- `nσ`: number of levels in the vertical
- `x0`: grid origin x co-ordinate 
- `y0`: grid origin y co-ordinate
- `h_mask`: Mask defining domain points within grid
- `h_isfixed': Mask defining locations of fixed thickness within grid
- `hyd_potential_isfixed': Mask defining locations of fixed hydraulic potential at the bed within grid
- `u_iszero`: Locations of zero u velocity points
- `v_iszero`: Locations of zero v velocity points
- `u_isfixed`: Locations of fixed u velocity points
- `v_isfixed`: Locations of fixed v velocity points
- `quadrature_weights`: weights associated with sigma levels used in quadrature scheme
- `basin_ID`: grid of basin IDs 
"""
function Grid(; 
    nx = 80,
    ny = 10,
    dx = 8000.0,
    dy = 8000.0,
    nσ = 4,
    x0 = 0.0,
    y0 = -40000.0,
    h_mask = nothing,
    h_isfixed = nothing,
    hyd_potential_isfixed = nothing,
    u_iszero = nothing,
    v_iszero = nothing,
    u_isfixed = nothing,
    v_isfixed = nothing,
    quadrature_weights = nothing,
    σ = nothing,
    basin_ID = nothing)

#check integer inputs
((typeof(nx) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in x direction (nx) must a positive integer larger than one")) 
((typeof(ny) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in y direction (ny)  must a positive integer larger than one")) 
((typeof(nσ) <: Integer) && nx > 1) || throw(ArgumentError("number of grid cells in vertical (nσ)  must a positive integer larger than one")) 

#if boundary conditions passed as string array, assemble these matric
~(typeof(u_iszero) == Vector{String}) || (u_iszero = orientations2bc(deepcopy(u_iszero),nx+1,ny))
~(typeof(v_iszero) == Vector{String}) || (v_iszero = orientations2bc(deepcopy(v_iszero),nx,ny+1))
~(typeof(h_isfixed) == Vector{String}) || (h_isfixed = orientations2bc(deepcopy(h_isfixed),nx,ny))
~(typeof(hyd_potential_isfixed) == Vector{String}) || (hyd_potential_isfixed = orientations2bc(deepcopy(hyd_potential_isfixed),nx,ny))
~(typeof(u_isfixed) == Vector{String}) || (u_isfixed = orientations2bc(deepcopy(u_isfixed),nx+1,ny))
~(typeof(v_isfixed) == Vector{String}) || (v_isfixed = orientations2bc(deepcopy(v_isfixed),nx,ny+1))


#assemble h_mask, u_iszero, v_iszero (if not passed as string)
(~(h_mask === nothing)) || (h_mask = trues(nx,ny))
(~(h_isfixed === nothing)) || (h_isfixed = falses(nx,ny))
(~(hyd_potential_isfixed === nothing)) || (hyd_potential_isfixed = falses(nx,ny))
(~(u_iszero === nothing))|| (u_iszero = falses(nx+1,ny))
(~(v_iszero === nothing)) || (v_iszero = falses(nx, ny+1))
(~(u_isfixed === nothing))|| (u_isfixed = falses(nx+1,ny))
(~(v_isfixed === nothing)) || (v_isfixed = falses(nx, ny+1))
(~(basin_ID === nothing)) || (basin_ID = ones(nx,ny))

#check the sizes of inputs
size(h_mask)==(nx,ny) || throw(DimensionMismatch("h_mask size must be (nx x ny) (i.e. $nx x $ny)"))
size(h_isfixed)==(nx,ny) || throw(DimensionMismatch("h_isfixed size must be (nx x ny) (i.e. $nx x $ny)"))
size(hyd_potential_isfixed)==(nx,ny) || throw(DimensionMismatch("hyd_potential_isfixed size must be (nx x ny) (i.e. $nx x $ny)"))
size(u_iszero)==(nx+1,ny) || throw(DimensionMismatch("u_iszero size must be size of U grid (nx+1 x ny) (i.e. $(nx+1) x $ny)"))
size(v_iszero)==(nx,ny+1) || throw(DimensionMismatch("v_iszero size must be size of V grid (nx x ny+1) (i.e. $nx x $(ny+1)"))
size(basin_ID)==(nx,ny) || throw(DimensionMismatch("Basin_ID size must be (nx x ny) (i.e. $nx x $ny)"))
size(u_isfixed)==(nx+1,ny) || throw(DimensionMismatch("u_isfixed size must be size of U grid (nx+1 x ny) (i.e. $(nx+1) x $ny)"))
size(v_isfixed)==(nx,ny+1) || throw(DimensionMismatch("v_isfixed size must be size of V grid (nx x ny+1) (i.e. $nx x $(ny+1)"))


#map bit arrays to boolean
try
    h_mask = convert(Array{Bool,2}, h_mask)
    #@assert h_mask == clip(h_mask)
catch 
    throw(ArgumentError("h_mask must be Boolean (or equivalent)"))
end

try
    h_isfixed = convert(Array{Bool,2}, h_isfixed)
catch 
    throw(ArgumentError("h_isfixed must be Boolean (or equivalent)"))
end

try
    hyd_potential_isfixed = convert(Array{Bool,2}, hyd_potential_isfixed)
catch 
    throw(ArgumentError("hyd_potential_isfixed must be Boolean (or equivalent)"))
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

try 
    u_isfixed = convert(Array{Bool,2}, u_isfixed)
catch 
    throw(ArgumentError("u_isfixed must be Boolean (or equivalent)"))
end

try
    v_isfixed = convert(Array{Bool,2}, v_isfixed)
catch
    throw(ArgumentError("v_isfixed must be Boolean (or equivalent)"))
end

#compute grid co-ordinates
xxh=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xxh)==(nx,ny)
yyh=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yyh)==(nx,ny)

xxu=[x0+(i-1.0)*dx for i=1:(nx+1), j=1:ny]; @assert size(xxu)==(nx+1,ny)
yyu=[y0+(j-0.5)*dy for i=1:(nx+1), j=1:ny]; @assert size(yyu)==(nx+1,ny)

xxv=[x0+(i-0.5)*dx for i=1:nx, j=1:(ny+1)]; @assert size(xxv)==(nx,ny+1)
yyv=[y0+(j-1.0)*dy for i=1:nx, j=1:(ny+1)]; @assert size(yyv)==(nx,ny+1)

xxc=[x0+i*dx for i=1:(nx-1), j=1:(ny-1)]; @assert size(xxc)==(nx-1,ny-1)
yyc=[y0+j*dy for i=1:(nx-1), j=1:(ny-1)]; @assert size(yyc)==(nx-1,ny-1)

#sigma grid info and checks
~(σ === nothing) ||  (σ = collect(range(0.0,length=nσ,stop=1.0))) #default sigma
length(σ) == nσ || throw(DimensionMismatch("number of sigma levels (= $(length(σ))) must match number of σ grid points (=$(nσ))"))
ζ = one(eltype(σ)) .- σ ; @assert length(ζ) == nσ

@assert σ[1] == zero(eltype(σ)) && σ[end] == one(eltype(σ))
(~(quadrature_weights === nothing) || (quadrature_weights = 0.5*[ σ[2] .- σ[1] ; σ[3:end] .- σ[1:end-2] ; σ[end] .- σ[end-1] ]))
size(quadrature_weights) == (nσ,) || throw(DimensionMismatch("Input quadrate weighs are size $size(quadrature_weights). quadrature weights must have size (nσ,) (i.e. ($nσ,))"))


return Grid(nx,
            ny,
            nσ,
            dx,
            dy,
            x0,
            y0,
            h_mask,
            h_isfixed,
            hyd_potential_isfixed,
            u_iszero,
            v_iszero,
            u_isfixed,
            v_isfixed,
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
            quadrature_weights,
            basin_ID)
end

"""
    orientations2bc(directions, M, N)

Make an M x N matrix with trues in the locations specified by directions
"""
function orientations2bc(orientations, M, N)
    A = falses(M,N)
    for ornt in orientations
        ornt_low = lowercase(ornt)
        ~(ornt_low == "north") || (A[1,:] .= true)
        ~(ornt_low == "south") || (A[end,:] .= true)
        ~(ornt_low == "east") || (A[:,end] .= true)
        ~(ornt_low == "west") || (A[:,1] .= true)

    end
    return A 
end

Base.size(g::Grid) = (g.nx, g.ny, g.nσ)

function Base.show(io::IO, g::Grid)
    return print(io, "Grid ", summary(g))
end


"""
    reconstruct_on_grid(s, grid::Grid)

Either return s unchanged, or return a new instance that is modified for use on grid.
Different types can overload this to implement their own specialised versions.

"""
function reconstruct_on_grid(s, grid::Grid)
    return s
end

"""
    reconstruct_on_grid(s, params::Params, grid::Grid)

Either return s unchanged, or return a new instance that is modified for use on grid.
Different types can overload this to implement their own specialised versions.

"""
function reconstruct_on_grid(s, params, grid::Grid)
    return s
end


"""
    spatial_on_subdomain(a, grid, subdomain)

Copy the rectangle of `a` that belongs to one piece of `grid`.

`subdomain` is `(x_start, x_end, y_start, y_end)`, the index box of one MPI rank or
Schwarz tile. Numbers, functions, and `nothing` are left unchanged. A 2D or 3D array
that matches the ice grid in x and y is copied on that box (3D keeps every vertical
level). Other sizes are left unchanged, so a small placeholder array is not sliced.
"""
function spatial_on_subdomain(a, grid::Grid, subdomain::NTuple{4,<:Integer})
    (isnothing(a) || a isa Number || a isa Function) && return a
    x_start, x_end, y_start, y_end = subdomain
    if ndims(a) == 2 && size(a) == (grid.nx, grid.ny)
        return copy(a[x_start:x_end, y_start:y_end])
    elseif ndims(a) == 3 && size(a, 1) == grid.nx && size(a, 2) == grid.ny
        return copy(a[x_start:x_end, y_start:y_end, axes(a, 3)])
    else
        return a
    end
end

"""
    field_on_grid(a, grid; name=nothing)

Turn a number into an array of that value on every cell of `grid`. Leave a matching
2D array as it is.

If you pass `name` and the array is the wrong size, throw `DimensionMismatch` with
that name in the message. Without `name`, a wrong-sized array is left unchanged.
"""
function field_on_grid(a, grid::Grid; name::Union{Nothing,AbstractString} = nothing)
    if a isa Number
        return a .* ones(grid.nx, grid.ny)
    elseif name === nothing || (ndims(a) == 2 && size(a) == (grid.nx, grid.ny))
        return a
    else
        throw(DimensionMismatch("$name does not match grid size $(grid.nx) x $(grid.ny)"))
    end
end

"""
    subdomain_index_ranges(x_indices, y_indices, grid, subdomain)

Cut the file-lookup lists down to one MPI rank or Schwarz tile of `grid`.

`x_indices` and `y_indices` say, for each local cell, which column to read from a
global file. If they are `nothing`, start from `1:nx` and `1:ny`. If they already
exist, they are cut to `subdomain` rather than started again from 1.
"""
function subdomain_index_ranges(x_indices, y_indices, grid::Grid, subdomain::NTuple{4,<:Integer})
    x_start, x_end, y_start, y_end = subdomain
    parent_x = isnothing(x_indices) ? (1:grid.nx) : x_indices
    parent_y = isnothing(y_indices) ? (1:grid.ny) : y_indices
    return parent_x[x_start:x_end], parent_y[y_start:y_end]
end

"""
    grid_index_ranges(x_indices, y_indices, grid)

File-lookup lists to use on `grid` (an MPI/Schwarz tile, or the full domain).

Keep lists if they already exist (this object is already a tile). Otherwise use
`1:nx` and `1:ny`, so local index 1 is file index 1 (the full-domain case).
"""
function grid_index_ranges(x_indices, y_indices, grid::Grid)
    xs = isnothing(x_indices) ? (1:grid.nx) : x_indices
    ys = isnothing(y_indices) ? (1:grid.ny) : y_indices
    return xs, ys
end

"""
    forcing_index_ranges(x_indices, y_indices, dest)

Which columns of a forcing file to read into local array `dest` (MPI tile or full domain).

Use the lookup lists if they exist. If they are `nothing`, read a `dest`-sized block
starting at 1. The list lengths must match `dest` in x and y.
"""
function forcing_index_ranges(x_indices, y_indices, dest)
    is = isnothing(x_indices) ? (1:size(dest, 1)) : x_indices
    js = isnothing(y_indices) ? (1:size(dest, 2)) : y_indices
    (length(is) == size(dest, 1) && length(js) == size(dest, 2)) || throw(DimensionMismatch(
        "x/y index ranges ($(length(is)),$(length(js))) do not match field size $(size(dest)[1:2])"
    ))
    return is, js
end

"""
    reconstruct_on_subdomain(s, grid::Grid, subdomain::NTuple{4,<: Integer})

Either return s unchanged, or return a new instance that is modified for use on a subdomain of a grid.
The subdomain (i_start,i_end,j_start,j_end) contains indices for the portion of the grid [i_start:i_end,j_start:j_end]. 
Different types can overload this to implement their own specialised versions.

"""
function reconstruct_on_subdomain(s, grid::Grid, subdomain::NTuple{4,<: Integer}) 
    return s
end


end
