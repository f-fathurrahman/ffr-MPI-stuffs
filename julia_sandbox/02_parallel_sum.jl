using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

N = 10^7 ÷ size
data = rand(N)

MPI.Barrier(comm)
t0 = MPI.Wtime()

local_sum = sum(data)
global_sum = MPI.Allreduce(local_sum, +, comm)

MPI.Barrier(comm)
t1 = MPI.Wtime()

elapsed = t1 - t0
parallel_time = MPI.Reduce(elapsed, max, 0, comm)

if rank == 0
    println("Parallel time = $parallel_time seconds")
end

MPI.Finalize()


