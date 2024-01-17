all:sparc

CXX = mpicxx.openmpi

PETSC_DIR1 = /home/efefer/WORKS/PETSC/petsc-3.5.3
#PETSC_DIR2 = /home/efefer/WORKS/PETSC/petsc-3.5.3/arch-linux2-c-debug
PETSC_DIR2 = /home/efefer/WORKS/PETSC/petsc-3.5.3/arch-linux2-c-opt

MKLROOT = /home/efefer/mysoftwares/LIBMKL/

CPPFLAGS = -O3 -I ${PETSC_DIR1}/include -I ${PETSC_DIR2}/include -I ./inc -I ${MKLROOT}/include

OBJSC = main.o initObjs.o spline.o nonlocal.o density.o \
readfiles.o energy.o ExchangeCorrelation.o occupation.o poisson.o \
chebyshev.o scf.o mixing.o forces.o relaxatoms.o multipole.o kPointHamiltonian.o

PETSC_LIB = -Wl,-rpath,${PETSC_DIR2}/lib: ${PETSC_DIR2}/lib/libpetsc.so \
-L${MKLROOT} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

CLEANFILES = ./lib/sparc ${OBJSC}

sparc: ${OBJSC}
	${CXX} -Wall -o sparc.x ${OBJSC} ${PETSC_LIB}

