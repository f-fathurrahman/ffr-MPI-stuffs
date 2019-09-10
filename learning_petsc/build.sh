PETSC_DIR1=/home/efefer/WORKS/PETSC/petsc-3.11.3/
PETSC_DIR2=/home/efefer/WORKS/PETSC/petsc-3.11.3/arch-linux2-c-debug

mpicc.mpich $1 -I${PETSC_DIR1}/include -I${PETSC_DIR2}/include -L${PETSC_DIR2}/lib -lpetsc

# -lblas -llapack -lm -lX11 -ldl
