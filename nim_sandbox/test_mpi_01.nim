{.passL: "-lmpi".}   # Link against MPI (optional if using mpicc)

type
  MPI_Comm* = distinct int32
  MPI_Datatype* = distinct int32
  MPI_Status* = object
    MPI_SOURCE*: int32
    MPI_TAG*: int32
    MPI_ERROR*: int32

const
  MPI_COMM_WORLD* {.importc.}: MPI_Comm
  MPI_INT* {.importc.}: MPI_Datatype

proc MPI_Init(argc: ptr cint; argv: ptr ptr cchar): cint {.importc.}
proc MPI_Finalize(): cint {.importc.}
proc MPI_Comm_size(comm: MPI_Comm; size: ptr cint): cint {.importc.}
proc MPI_Comm_rank(comm: MPI_Comm; rank: ptr cint): cint {.importc.}
proc MPI_Send(buf: pointer; count: cint; datatype: MPI_Datatype;
              dest: cint; tag: cint; comm: MPI_Comm): cint {.importc.}
proc MPI_Recv(buf: pointer; count: cint; datatype: MPI_Datatype;
              source: cint; tag: cint; comm: MPI_Comm;
              status: ptr MPI_Status): cint {.importc.}

when isMainModule:
  var argc = 0.cint
  var argv: ptr ptr cchar = nil

  discard MPI_Init(addr argc, addr argv)

  var size, rank: cint
  discard MPI_Comm_size(MPI_COMM_WORLD, addr size)
  discard MPI_Comm_rank(MPI_COMM_WORLD, addr rank)

  echo "Hello from rank ", rank, " of ", size

  discard MPI_Finalize()
  

