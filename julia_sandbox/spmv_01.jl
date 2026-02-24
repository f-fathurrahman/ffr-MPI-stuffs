using MPI

function main()

    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    ############################
    # GLOBAL GRID
    ############################

    Nx, Ny, Nz = 64, 64, 64

    ############################
    # PROCESS GRID
    ############################

    #dims = MPI.Dims_create(nprocs, 3)
    #dims = [0,0,0]                  # must be a Vector{Int}
    #MPI.Dims_create!(nprocs, 3, dims)    # fills dims in-place
    #cart = MPI.Cart_create(comm, dims; reorder=true, periods=(false,false,false))
    #periods = [0, 0, 0]
    #reorder = true
    #cart = MPI.Cart_create(comm, dims, periods, reorder )

    ndims = 3
    dims = MPI.Dims_create(nprocs, zeros(Int64, ndims))

    periods = ntuple(_ -> false, ndims)
    cart = MPI.Cart_create(comm, dims, periodic=periods, reorder=true)

    coords = MPI.Cart_coords(cart, rank)

    Px, Py, Pz = dims
    ix, iy, iz = coords

    nx = div(Nx, Px)
    ny = div(Ny, Py)
    nz = div(Nz, Pz)

    ############################
    # LOCAL ARRAY WITH HALO
    ############################

    u = zeros(Float64, nx+2, ny+2, nz+2)
    y = zeros(Float64, nx, ny, nz)

    # fill interior with test data
    for k in 1:nz, j in 1:ny, i in 1:nx
        u[i+1,j+1,k+1] = rank + 1
    end

    ############################
    # NEIGHBORS
    ############################

    function neighbor(dir, disp)
        src, dst = MPI.Cart_shift(cart, dir, disp)
        return dst
    end

    nbr_xm = neighbor(1, -1)
    nbr_xp = neighbor(1, +1)
    nbr_ym = neighbor(2, -1)
    nbr_yp = neighbor(2, +1)
    nbr_zm = neighbor(3, -1)
    nbr_zp = neighbor(3, +1)

    ############################
    # HALO EXCHANGE
    ############################

    function exchange_halo!(u)

        reqs = MPI.Request[]

        ################ X faces ################
        sendbuf = view(u,2,:,:)
        recvbuf = view(u,1,:,:)
        if nbr_xm != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_xm, 0, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_xm, 1, cart))
        end

        sendbuf = view(u,nx+1,:,:)
        recvbuf = view(u,nx+2,:,:)
        if nbr_xp != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_xp, 1, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_xp, 0, cart))
        end

        ################ Y faces ################
        sendbuf = view(u,:,2,:)
        recvbuf = view(u,:,1,:)
        if nbr_ym != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_ym, 2, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_ym, 3, cart))
        end

        sendbuf = view(u,:,ny+1,:)
        recvbuf = view(u,:,ny+2,:)
        if nbr_yp != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_yp, 3, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_yp, 2, cart))
        end

        ################ Z faces ################
        sendbuf = view(u,:,:,2)
        recvbuf = view(u,:,:,1)
        if nbr_zm != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_zm, 4, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_zm, 5, cart))
        end

        sendbuf = view(u,:,:,nz+1)
        recvbuf = view(u,:,:,nz+2)
        if nbr_zp != MPI.PROC_NULL
            push!(reqs, MPI.Isend(sendbuf, nbr_zp, 5, cart))
            push!(reqs, MPI.Irecv!(recvbuf, nbr_zp, 4, cart))
        end

        MPI.Waitall(reqs)
    end

    ############################
    # LAPLACIAN APPLY
    ############################

    function apply_laplacian!(y,u)

        for k in 1:nz, j in 1:ny, i in 1:nx
            ii = i+1; jj=j+1; kk=k+1

            y[i,j,k] =
                6*u[ii,jj,kk] -
                u[ii-1,jj,kk] -
                u[ii+1,jj,kk] -
                u[ii,jj-1,kk] -
                u[ii,jj+1,kk] -
                u[ii,jj,kk-1] -
                u[ii,jj,kk+1]
        end
    end

    ############################
    # RUN
    ############################

    exchange_halo!(u)
    apply_laplacian!(y,u)

    ############################
    # CHECKSUM
    ############################

    local_sum = sum(y)
    global_sum = MPI.Allreduce(local_sum, +, cart)

    if rank==0
        println("Global checksum = ", global_sum)
    end

    MPI.Finalize()

end

main()

