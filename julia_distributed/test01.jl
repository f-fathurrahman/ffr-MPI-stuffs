using Distributed

# Add worker processes (adjust as needed)
addprocs(4)   # or Sys.CPU_THREADS - 1

@everywhere using Random

@info nprocs()
@info workers()
