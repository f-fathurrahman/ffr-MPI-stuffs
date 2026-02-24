using MPI

# works for nprocs >= 3

# ------------------------------------------------------------
# Uneven domain decomposition (like PETSc DMDA)
# ------------------------------------------------------------
function local_size(N, P, coord)
    base = div(N, P)
    remainder = N % P
    return coord < remainder ? base + 1 : base
end

function global_offset(N, P, coord)
    base = div(N, P)
    remainder = N % P

    if coord < remainder
        return coord * (base + 1)
    else
        return remainder * (base + 1) +
               (coord - remainder) * base
    end
end


function halo_exchange!(u, cart_comm, dims)

    ndims = 3

    for dim in 0:ndims-1

        # Skip if dimension not distributed
        if dims[dim+1] == 1
            continue
        end

        src, dst = MPI.Cart_shift(cart_comm, dim, 1)

        if dim == 0
            MPI.Sendrecv!(view(u, 2,:,:), dst, 10,
                          view(u, 1,:,:), src, 10,
                          cart_comm)

            MPI.Sendrecv!(view(u, size(u,1)-1,:,:), src, 11,
                          view(u, size(u,1),:,:), dst, 11,
                          cart_comm)

        elseif dim == 1
            MPI.Sendrecv!(view(u, :,2,:), dst, 20,
                          view(u, :,1,:), src, 20,
                          cart_comm)

            MPI.Sendrecv!(view(u, :,size(u,2)-1,:), src, 21,
                          view(u, :,size(u,2),:), dst, 21,
                          cart_comm)

        elseif dim == 2
            MPI.Sendrecv!(view(u, :,:,2), dst, 30,
                          view(u, :,:,1), src, 30,
                          cart_comm)

            MPI.Sendrecv!(view(u, :,:,size(u,3)-1), src, 31,
                          view(u, :,:,size(u,3)), dst, 31,
                          cart_comm)
        end
    end
end


# ------------------------------------------------------------
# 7-point Laplacian (SPMV-like stencil)
# ------------------------------------------------------------
function laplacian!(Lu, u)
    nx, ny, nz = size(u)

    for i in 2:nx-1, j in 2:ny-1, k in 2:nz-1
        Lu[i,j,k] =
            -6.0*u[i,j,k] +
             u[i+1,j,k] + u[i-1,j,k] +
             u[i,j+1,k] + u[i,j-1,k] +
             u[i,j,k+1] + u[i,j,k-1]
    end
end


# ------------------------------------------------------------
# Main program
# ------------------------------------------------------------
function main()

    MPI.Init()

    world = MPI.COMM_WORLD
    rank = MPI.Comm_rank(world)
    nprocs = MPI.Comm_size(world)

    # --------------------------------------------------------
    # Create Cartesian topology (works for ANY nprocs)
    # --------------------------------------------------------
    ndims = min(3, nprocs)
    dims = MPI.Dims_create(nprocs, zeros(Int64,ndims))

    periods = ntuple(_ -> false, ndims)

    cart_comm = MPI.Cart_create(world, dims, periodic=periods, reorder=true)
    coords = MPI.Cart_coords(cart_comm, rank)

    if rank == 0
        println("Processor grid = $dims")
    end

    # --------------------------------------------------------
    # Global grid size
    # --------------------------------------------------------
    Nx, Ny, Nz = 64, 64, 64

    # Handle lower dimensional cases automatically
    dims_full = ntuple(i -> i <= ndims ? dims[i] : 1, 3)
    coords_full = ntuple(i -> i <= ndims ? coords[i] : 0, 3)

    nx = local_size(Nx, dims_full[1], coords_full[1])
    ny = local_size(Ny, dims_full[2], coords_full[2])
    nz = local_size(Nz, dims_full[3], coords_full[3])

    # --------------------------------------------------------
    # Allocate local grid with ghost layers
    # --------------------------------------------------------
    u  = rand(nx+2, ny+2, nz+2)
    Lu = zeros(nx+2, ny+2, nz+2)

    println("Rank $rank coords=$coords_full local=($nx,$ny,$nz)")

    # --------------------------------------------------------
    # Halo exchange
    # --------------------------------------------------------
    halo_exchange!(u, cart_comm, dims)

    # --------------------------------------------------------
    # Apply Laplacian (SPMV)
    # --------------------------------------------------------
    laplacian!(Lu, u)

    local_norm = sum(abs2, Lu)
    global_norm = MPI.Allreduce(local_norm, +, world)

    if rank == 0
        println("Global norm = $global_norm")
    end

    MPI.Barrier(world)
    MPI.Finalize()
end


main()
