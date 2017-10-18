# MPI functions

The old way:

```fortran
PROGRAM my_mpi_program
  IMPLICIT NONE
  include 'mpif.h'
  ! rest of the program  ...
END PROGRAM
```

Initialize MPI

```fortran
INTEGER :: ierr
CALL MPI_Init( ierr )
```

Finalize MPI

```fortran
INTEGER :: ierr
CALL MPI_Finalize( ierr )
```

Get processor name

```fortran
INTEGER :: ierr
CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname
INTEGER :: length
CALL MPI_Get_processor_name( hostname, length, ierr )
```

