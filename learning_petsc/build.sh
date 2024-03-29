#!/bin/bash

PETSC_DIR=/opt/petsc-3.13.0_openmpi-2.1.1_debug

LINK1="-Wl,-rpath,$PETSC_DIR/lib: $PETSC_DIR/lib/libpetsc.so"

if [ "$#" -eq 1 ]; then

  basnam=`basename $1 .c`
  #mpicc.openmpi $1 -I${PETSC_DIR}/include -L${PETSC_DIR}/lib -lpetsc -o $basnam.x -lm
  mpicc.openmpi $1 -I${PETSC_DIR}/include -o $basnam.x $LINK1 -lm
  # -lblas -llapack -lm -lX11 -ldl

else

  echo "Wrong number of arguments: $#"

fi

