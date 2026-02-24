using MPI

function main()

    MPI.Init()

    world = MPI.COMM_WORLD
    rank = MPI.Comm_rank(world)
    nprocs = MPI.Comm_size(world)

    ndims = 3
    dims = MPI.Dims_create(nprocs, zeros(Int64, ndims))

    periods = ntuple(_ -> false, ndims)
    cart_comm = MPI.Cart_create(world, dims, periodic=periods, reorder=true)

    coords = MPI.Cart_coords(cart_comm, rank)

    if rank == 0
        println("Processor grid = $dims")
    end

    println("Rank $rank coords = $coords")

    MPI.Finalize()
end

main()