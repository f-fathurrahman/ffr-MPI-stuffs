! adapted from LLNL MPI exercises
PROGRAM hello

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ierr, taskid, numtasks, length
  CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname
  
  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, numtasks, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, taskid, ierr )

  CALL MPI_Get_processor_name( hostname, length, ierr )

  WRITE(*,*) 'Hello from: taskid = ', taskid, ' hostname ', trim(hostname)
  IF( taskid == 0 ) THEN
    WRITE(*,*) 'MPI_MAX_PROCESSOR_NAME = ', MPI_MAX_PROCESSOR_NAME
    WRITE(*,*) 'Number of tasks = ', numtasks
  ENDIF
 
  CALL MPI_Finalize( ierr )

END PROGRAM 
