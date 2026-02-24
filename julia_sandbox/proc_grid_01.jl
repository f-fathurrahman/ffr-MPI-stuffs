using MPI

function main()

    MPI.Init()
    comm = MPI.COMM_WORLD
    Nprocs = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    # Define k-point and band parallelization
    Nk = 2   # number of k-point groups
    Nb = Nprocs ÷ Nk  # bands per k

    if rank == 0
        println("Nprocs = $Nprocs")
        println("Nk = $Nk and Nb = $Nb")
    end

    @assert Nk * Nb == Nprocs

    # -----------------------------------
    # Split by k-point
    # -----------------------------------

    k_color = rank ÷ Nb
    k_comm = MPI.Comm_split(comm, k_color, rank)

    k_rank = MPI.Comm_rank(k_comm)
    k_size = MPI.Comm_size(k_comm)

    # -----------------------------------
    # Split inside k_comm by band
    # -----------------------------------

    band_color = k_rank
    band_comm = MPI.Comm_split(comm, band_color, rank)

    band_rank = MPI.Comm_rank(band_comm)
    band_size = MPI.Comm_size(band_comm)

    println("Global rank=$rank | k_group=$k_color (k_rank=$k_rank) | band_rank=$band_rank")

    # -----------------------------------
    # Dummy workload: wavefunction norms
    # -----------------------------------

    local_wavefunc = rand(10_000) .* (rank + 1)

    local_norm = sum(abs2, local_wavefunc)

    # Reduce within band communicator (like orthonormalization)
    band_norm = MPI.Allreduce(local_norm, +, band_comm)

    # Reduce within k communicator (like energy per k-point)
    k_energy = MPI.Allreduce(local_norm, +, k_comm)

    if k_rank == 0
        println("k_group $k_color total energy = $k_energy")
    end

    MPI.Barrier(comm)
    MPI.Finalize()
end

main()

