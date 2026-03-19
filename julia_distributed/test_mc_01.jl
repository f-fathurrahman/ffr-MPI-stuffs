using Distributed

# Add worker processes (adjust as needed)
addprocs(4)   # or Sys.CPU_THREADS - 1

@everywhere using Random

@info nprocs()
@info workers()

function mc_pi_serial(N)
    count = 0
    for i in 1:N
        x = rand()
        y = rand()
        count += (x^2 + y^2 <= 1)
    end
    return 4 * count / N
end

@everywhere function mc_pi_chunk(N)
    count = 0
    for i in 1:N
        x = rand()
        y = rand()
        count += (x^2 + y^2 <= 1)
    end
    return count
end

function mc_pi_parallel(N_total)
    nprocs = nworkers()
    N_per_worker = div(N_total, nprocs)

    counts = pmap(mc_pi_chunk, fill(N_per_worker, nprocs))
    total_count = sum(counts)

    return 4 * total_count / (N_per_worker * nprocs)
end


N = 10^8
@time mc_pi_serial(N)
@time mc_pi_parallel(N)

