using MPI
using WAVI
using WAVI.Specs

export @root

"""
    Get directional halos for a given MPISpec

    This function returns the halos for a given MPISpec in the form of a tuple (top, right, bottom, left).
    If a neighbour is not present, the halo is set to 0.

    Args:
        spec: The MPISpec to get the halos for

    Returns:
        A tuple of the form (top, right, bottom, left)
"""
function get_halos(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    # Note: MPI_Cart_shift returns `MPI_PROC_NULL` (represented as a negative integer)
    #       if the neighbour is not present
    return spec.top > -1 ? spec.halo : 0, 
           spec.right > -1 ? spec.halo : 0, 
           spec.bottom > -1 ? spec.halo : 0, 
           spec.left > -1 ? spec.halo : 0
end

"""
True when this rank has a cardinal MPI neighbour.

One-rank runs have no neighbours, so halo, PoU, and Schwarz stay uncompiled.
"""
mpi_has_neighbours(spec::MPISpec) =
    spec.left > -1 || spec.right > -1 || spec.top > -1 || spec.bottom > -1

"""
    Get the bounds of a given MPISpec

    The bounds are calculated based on the global grid size and the number of processes in each direction.

    Args:
        spec: The MPISpec to get the bounds for

    Returns:
        A tuple of the form (x_start, x_end, y_start, y_end)
"""
function get_bounds(spec::MPISpec)::Tuple{Int, Int, Int, Int}
    # Calculate base local grid size assuming even distribution
    nx_base, nx_rem = divrem(spec.global_grid.nx, spec.px)
    ny_base, ny_rem = divrem(spec.global_grid.ny, spec.py)

    # Adding extra cell to local grid size if there is a remainder
    local_nx = nx_base + (spec.coords[1] < nx_rem ? 1 : 0)
    local_ny = ny_base + (spec.coords[2] < ny_rem ? 1 : 0)

    # Calculate start index (accounting for previous ranks that took an extra cell)
    # (Rank * BaseSize) + (How many previous ranks took an extra cell) + 1
    x_start = spec.coords[1] * nx_base + min(spec.coords[1], nx_rem) + 1
    y_start = spec.coords[2] * ny_base + min(spec.coords[2], ny_rem) + 1

    # Calculate end index
    x_end = x_start + local_nx - 1
    y_end = y_start + local_ny - 1

    halos = get_halos(spec)
    return max(x_start - halos[4], 1),                   # Left: Protect against indices < 1
           min(x_end   + halos[2], spec.global_grid.nx), # Right: Protect against indices > global max
           max(y_start - halos[1], 1),                   # Bottom: Protect against indices < 1
           min(y_end   + halos[3], spec.global_grid.ny)  # Top: Protect against indices > global max
end


"""
    Get the size of a given MPISpec

    The size is calculated based on the global grid size and the number of processes in each direction.

    Args:
        spec: The MPISpec to get the size for

    Returns:
        A tuple of the form (nx, ny)
"""
# Get size of local grid
function get_size(spec::MPISpec)::Tuple{Int, Int}
    nx_base, nx_rem = divrem(spec.global_grid.nx, spec.px)
    ny_base, ny_rem = divrem(spec.global_grid.ny, spec.py)

    # If this rank's x-coord is < remainder, it gets an extra cell
    # This accounts for the fact that the global grid may not be evenly split
    # (e.g. if the global grid is 10x10 and we have 3 ranks in x, then the first
    # rank will have 4 cells in x and the other two ranks will have 3 cells in x)
    local_nx = nx_base + (spec.coords[1] < nx_rem ? 1 : 0)
    local_ny = ny_base + (spec.coords[2] < ny_rem ? 1 : 0)

    halos = get_halos(spec)
    local_nx += halos[4] + halos[2]
    local_ny += halos[1] + halos[3]

    return local_nx, local_ny
end

# Sourced from Oceananigans.DistributedComputations
# - available under the MIT License: https://github.com/CliMA/Oceananigans.jl/blob/main/LICENSE
# - original copyright of the below code (except adaptations for use) 2018 Climate Modeling Alliance

mpi_initialized()     = MPI.Initialized()
mpi_rank(comm)        = MPI.Comm_rank(comm)
mpi_size(comm)        = MPI.Comm_size(comm)
global_barrier(comm)  = MPI.Barrier(comm)
global_communicator() = MPI.COMM_WORLD

"""
    @root communicator exp...

Perform `exp` only on rank 0 in communicator, otherwise known as the "root" rank.
Other ranks will wait for the root rank to finish before continuing.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro root(communicator, exp)
    command = quote
        if WAVI.Specs.mpi_initialized()
            rank = WAVI.Specs.mpi_rank($communicator)
            if rank == 0
                $exp
            end
            WAVI.Specs.global_barrier($communicator)
        else
            $exp
        end
    end
    return esc(command)
end

macro root(exp)
    command = quote
        @root WAVI.Specs.global_communicator() $exp
    end
    return esc(command)
end

## END Sourced from Oceananigans.DistributedComputations