using Printf
import MPI

function main()
    MPI.Init()
    
    Nprocs = MPI.Comm_size(MPI.COMM_WORLD)
    my_rank = MPI.Comm_rank(MPI.COMM_WORLD)
    
    @printf("Hello, I am proc: %d from %d total processes\n", my_rank, Nprocs)

    MPI.Finalize()
end

main()
