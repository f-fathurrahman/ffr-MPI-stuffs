#!/bin/bash
#echo "Number of arguments: $#"

if [ "$#" -eq 2 ]; then

  export OMPI_FC=$1
  basnam=`basename $2 .f90`
  mpifort $2 -o $basnam.x

elif [ "$#" -eq 1 ]; then

  basnam=`basename $1 .f90`
  mpifort $1 -o $basnam.x

else

  echo "Wrong number of arguments: $#"

fi


