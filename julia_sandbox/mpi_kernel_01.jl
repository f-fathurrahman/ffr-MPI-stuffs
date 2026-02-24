using MPI

function main()

    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    outer = 20_000_000 ÷ size
    inner = 50

    MPI.Barrier(comm)
    t0 = MPI.Wtime()

    s = 0.0
    for i in 1:outer
        x = i + rank
        for j in 1:inner
            x = sin(x) + cos(x) + sqrt(abs(x))
        end
        s += x
    end

    total = MPI.Reduce(s, +, 0, comm)

    MPI.Barrier(comm)
    t1 = MPI.Wtime()

    if rank == 0
        println("Processes: ", size)
        println("Time: ", t1 - t0)
    end

    MPI.Finalize()
end

main()
