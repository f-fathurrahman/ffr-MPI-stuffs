type
  MPI_Comm {.importc: "MPI_Comm", header: "mpi.h".} = object

var
  MPI_COMM_WORLD {.importc: "MPI_COMM_WORLD", header: "mpi.h".}: MPI_Comm

proc MPI_Init(argc: ptr cint; argv: ptr ptr ptr cchar): cint
  {.importc: "MPI_Init", header: "mpi.h".}

proc MPI_Finalize(): cint
  {.importc: "MPI_Finalize", header: "mpi.h".}

proc MPI_Comm_size(comm: MPI_Comm; size: ptr cint): cint
  {.importc: "MPI_Comm_size", header: "mpi.h".}

proc MPI_Comm_rank(comm: MPI_Comm; rank: ptr cint): cint
  {.importc: "MPI_Comm_rank", header: "mpi.h".}

when isMainModule:
  var argc = 0.cint
  var argv: ptr ptr cchar = nil

  discard MPI_Init(addr argc, addr argv)

  var size, rank: cint
  discard MPI_Comm_size(MPI_COMM_WORLD, addr size)
  discard MPI_Comm_rank(MPI_COMM_WORLD, addr rank)

  echo "Rank ", rank, " of ", size

  discard MPI_Finalize()
  
  

