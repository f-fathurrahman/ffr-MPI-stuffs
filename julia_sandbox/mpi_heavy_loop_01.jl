using MPI

function main()

    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    N = 200_000_000
    chunk = div(N, size)

    start = rank*chunk + 1
    stop  = (rank+1)*chunk

    MPI.Barrier(comm)
    t0 = MPI.Wtime()

    s = 0.0
    for i in start:stop
        x = sin(i)*cos(i)*sqrt(i)
        s += x
    end

    total = MPI.Reduce(s, +, 0, comm)

    MPI.Barrier(comm)
    t1 = MPI.Wtime()

    if rank == 0
        println("Processes: ", size)
        println("Result: ", total)
        println("Time: ", t1 - t0)
    end

    MPI.Finalize()

end

main()

