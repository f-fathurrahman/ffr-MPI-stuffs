#!/bin/bash
#echo "Number of arguments: $#"

if [ "$#" -eq 2 ]; then

  export OMPI_CC=$1
  basnam=`basename $2 .c`
  mpicc $2 -o $basnam.x

elif [ "$#" -eq 1 ]; then

  basnam=`basename $1 .c`
  mpicc $1 -o $basnam.x

else

  echo "Wrong number of arguments: $#"

fi


