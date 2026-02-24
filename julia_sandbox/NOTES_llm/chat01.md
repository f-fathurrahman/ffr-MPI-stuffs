Can you give some examples using MPI.jl, also showing quite large workload (but trivial) such that the parallel speedup can be measured

---

Below are **practical MPI.jl examples** that are:

* trivial mathematically
* computationally heavy
* easy to scale
* good for benchmarking speedup

All examples use **wall-clock timing** and can be run with different process counts.

---

# 0. Setup

Install MPI.jl and an MPI implementation:

```bash
julia -e 'using Pkg; Pkg.add("MPI")'
```

Run programs like:

```bash
mpiexec -n 1 julia script.jl
mpiexec -n 2 julia script.jl
mpiexec -n 4 julia script.jl
mpiexec -n 8 julia script.jl
```

---

# 1. Monte-Carlo π (Classic Parallel Benchmark)

Embarrassingly parallel → ideal scaling test.

### Script: `mpi_pi.jl`

```julia
using MPI
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
```

### Why this is good

* no communication except one reduction
* perfect for measuring raw parallel efficiency

---

# 2. Large Vector Compute (CPU-Bound)

Heavy floating-point workload.

### Script: `mpi_heavyloop.jl`

```julia
using MPI
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
```

---

# 3. Artificial Heavy Work Kernel (Customizable Load)

You can scale cost arbitrarily by increasing inner loop.

### Script: `mpi_kernel.jl`

```julia
using MPI
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
```

Increase `inner` to simulate heavier kernels.

---

# 4. Measuring Speedup

Run with multiple process counts and record times:

| processes | time |
| --------- | ---- |
| 1         | T₁   |
| 2         | T₂   |
| 4         | T₄   |
| 8         | T₈   |

Compute:

[
\text{Speedup}(p) = T_1 / T_p
]

Efficiency:

[
\text{Efficiency} = \frac{T_1}{p T_p}
]

---

# 5. Optional: Automated Benchmark Script

You can run:

```bash
for p in 1 2 4 8 16
do
   mpiexec -n $p julia mpi_pi.jl
done
```

---

# 6. Tips for Reliable Benchmarks

To get meaningful scaling:

✔ workload per process must be large
✔ avoid printing from all ranks
✔ measure time after MPI.Barrier
✔ run multiple times and average

---

# 7. What Each Example Tests

| Example     | Tests                     |
| ----------- | ------------------------- |
| Monte Carlo | pure parallel scaling     |
| Heavy loop  | floating-point throughput |
| Kernel      | compute-intensive scaling |

---

✅ **Best starting test:** Monte-Carlo π
→ simplest and shows near-ideal scaling if system is configured correctly.

---

If you want, I can also provide:

* a **strong vs weak scaling benchmark suite**
* a **communication-heavy MPI test**
* a **hybrid MPI + Threads.jl benchmark**

Just tell me which direction you want to explore.

---

