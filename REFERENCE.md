# MPI functions

The old way:

```fortran
PROGRAM my_mpi_program
  IMPLICIT NONE
  INCLUDE 'mpif.h'
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
MPI_COMM :: comm
CALL MPI_Comm_rank( comm, rank, ierr )
```

## Point-to-point communication

### Sending data

```fortran
<type> :: buf(*)
INTEGER :: count, datatype, dest, tag, comm, ierr
MPI_Send( buf, count, datatype, dest, tag, comm, ierr )
```

### Receiving data

```fortran
<type> :: buf(*)
INTEGER :: count, datatype, source, tag, comm, ierr
INTEGER :: stat(MPI_STATUS_SIZE)
MPI_Recv( buf, count, datatype, source, tag, comm, stat, ierr )
```

`stat` object can be used to determine parameters that have not been
fixed by the call's argument.

