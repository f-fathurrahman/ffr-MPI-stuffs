using MPI

function main()

    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    N_total = 2_000_000_000       # increase for heavier workload
    N_local = div(N_total, size)

    MPI.Barrier(comm)
    t0 = MPI.Wtime()

    count = 0
    for i in 1:N_local
        x = rand()
        y = rand()
        if x*x + y*y ≤ 1.0
            count += 1
        end
    end

    global_count = MPI.Reduce(count, +, 0, comm)

    MPI.Barrier(comm)
    t1 = MPI.Wtime()

    if rank == 0
        π_est = 4 * global_count / N_total
        println("Processes: ", size)
        println("π estimate: ", π_est)
        println("Time: ", t1 - t0, " seconds")
    end

    MPI.Finalize()

end

main()