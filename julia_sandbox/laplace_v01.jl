using MPI

function main()

    MPI.Init()
    world = MPI.COMM_WORLD
    rank = MPI.Comm_rank(world)
    nprocs = MPI.Comm_size(world)

    # ---------------------------
    # 3D processor grid
    # ---------------------------

    dims = MPI.Dims_create(nprocs, 3)
    periods = (false, false, false)

    cart_comm = MPI.Cart_create(world, dims, periods, reorder=false)
    coords = MPI.Cart_coords(cart_comm, rank)

    px, py, pz = dims
    ix, iy, iz = coords

    println("Rank $rank at coords ($ix,$iy,$iz)")

    # ---------------------------
    # Local grid size
    # ---------------------------

    Nx, Ny, Nz = 64, 64, 64

    nx = Nx ÷ px
    ny = Ny ÷ py
    nz = Nz ÷ pz

    # Add halo layers
    u = rand(nx+2, ny+2, nz+2)

    # ---------------------------
    # Halo exchange
    # ---------------------------

    function halo_exchange!(u, cart_comm)
        for dim in 0:2
            src, dst = MPI.Cart_shift(cart_comm, dim, 1)

            if src != MPI.PROC_NULL
                sendbuf = view(u, 2,:,:)
                recvbuf = view(u, 1,:,:)
                MPI.Sendrecv!(sendbuf, dst, 0,
                            recvbuf, src, 0,
                            cart_comm)
            end

            if dst != MPI.PROC_NULL
                sendbuf = view(u, size(u,1)-1,:,:)
                recvbuf = view(u, size(u,1),:,:)
                MPI.Sendrecv!(sendbuf, src, 1,
                            recvbuf, dst, 1,
                            cart_comm)
            end
        end
    end

    halo_exchange!(u, cart_comm)

    # ---------------------------
    # Apply 7-point Laplacian
    # ---------------------------

    function laplacian(u)
        nx, ny, nz = size(u)
        Lu = zeros(nx, ny, nz)

        for i in 2:nx-1, j in 2:ny-1, k in 2:nz-1
            Lu[i,j,k] =
                -6u[i,j,k] +
                u[i+1,j,k] + u[i-1,j,k] +
                u[i,j+1,k] + u[i,j-1,k] +
                u[i,j,k+1] + u[i,j,k-1]
        end
        return Lu
    end

    Lu = laplacian(u)

    local_norm = sum(abs2, Lu)
    global_norm = MPI.Allreduce(local_norm, +, world)

    if rank == 0
        println("Global Laplacian norm = $global_norm")
    end

    MPI.Barrier(world)
    MPI.Finalize()
end

main()
