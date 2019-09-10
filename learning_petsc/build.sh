#!/bin/bash
PETSC_DIR1=/home/efefer/WORKS/PETSC/petsc-3.11.3/
PETSC_DIR2=/home/efefer/WORKS/PETSC/petsc-3.11.3/arch-linux2-c-debug

if [ "$#" -eq 1 ]; then

  basnam=`basename $1 .c`
  mpicc.mpich $1 -I${PETSC_DIR1}/include -I${PETSC_DIR2}/include -L${PETSC_DIR2}/lib -lpetsc -o $basnam.x
  # -lblas -llapack -lm -lX11 -ldl

else

  echo "Wrong number of arguments: $#"

fi

