PROGRAM hello_Bsend

  IMPLICIT NONE 
  INCLUDE 'mpif.h'

  INTEGER :: ierr, numtasks, taskid, length
  INTEGER :: partner, message, stat(MPI_STATUS_SIZE)
  CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname
  INTEGER :: rem

  CALL MPI_Init( ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, taskid, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, numtasks, ierr )

  rem = mod( numtasks, 2 )

  IF ( rem /= 0 ) THEN 

    IF ( taskid == 0 ) THEN 
      WRITE(*,*) 'Number of tasks must be even: ', numtasks
    ENDIF

  ELSE

    CALL MPI_Get_processor_name( hostname, length, ierr )
    WRITE(*,*) 'taskid = ', taskid, ' hostname = ', trim(hostname)
    IF ( taskid == 0 ) THEN 
      WRITE(*,*) 'numtasks = ', numtasks
    ENDIF

    ! Determine partner and send/receive with partner
    IF ( taskid < numtasks/2 ) THEN 
      partner = numtasks/2 + taskid
      CALL MPI_Send( taskid, 1, MPI_INTEGER, partner, 1, MPI_COMM_WORLD, ierr )
      CALL MPI_Recv( message, 1, MPI_INTEGER, partner, 1, MPI_COMM_WORLD, stat, ierr )
    ELSEIF ( taskid >= numtasks/2 ) THEN
      partner = taskid - numtasks/2
      CALL MPI_Recv( message, 1, MPI_INTEGER, partner, 1, MPI_COMM_WORLD, stat, ierr )
      CALL MPI_Send( taskid, 1, MPI_INTEGER, partner, 1, MPI_COMM_WORLD, ierr )
    ENDIF

    WRITE(*,*) 'Task: ', taskid, ' is partner with: ', message 

  ENDIF  ! check for even number of procs

  CALL MPI_Finalize( ierr )

END PROGRAM 

