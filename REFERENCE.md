# MPI functions

The old way:

```fortran
PROGRAM my_mpi_program
  IMPLICIT NONE
  include 'mpif.h'
  ! rest of the program  ...
END PROGRAM
```

## Initialize MPI

```fortran
INTEGER :: ierr
CALL MPI_Init( ierr )
```

## Finalize MPI

```fortran
INTEGER :: ierr
CALL MPI_Finalize( ierr )
```

## Get processor name

```fortran
INTEGER :: ierr
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname
INTEGER :: length
CALL MPI_Get_processor_name( hostname, length, ierr )
```

## Getting total number of processors

```fortran
INTEGER :: ierr
INTEGER :: Nprocs
MPI_COMM :: comm  ! usually MPI_COMM_WORLD
CALL MPI_Comm_size( comm, Nprocs, ierr )
```

## Getting process' rank

```fortran
INTEGER :: ierr
INTEGER :: rank
MPI_COMM :: MPI_COMM_WORLD
CALL MPI_Comm_rank( comm, rank, ierr )
```

## Point-to-point communication

Sending data

```fortran
MPI_Send()
```

Receiving data

```fortran
MPI_Recv()
```

