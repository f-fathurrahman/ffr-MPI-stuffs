USE_MKL       = 0
USE_SCALAPACK = 0

# (USE_DP_SUBEIG = 1 is required if both USE_MKL = 0 and USE_SCALAPACK = 0)
USE_DP_SUBEIG = 1

# Set USE_FFTW = 1 to use FFTW for fast Fourier transform in vdWDF.
# Don't open it together with USE_MKL
USE_FFTW      = 1

DEBUG_MODE    = 1

ENABLE_SIMD_COMPLEX = 0

CPPFLAGS = -I./ -Ixc/exx/include -Ixc/vdW/d3/include -Ixc/vdW/vdWDF/include -Ixc/mgga/include -IhighT/include
LDLIBS   = -lrt
LDLIBS += -llapacke -llapack -lblas -lm

CPPFLAGS += -DUSE_FFTW
LDLIBS += -lfftw3_mpi -lfftw3 # fftw needs to be loaded before compile

# To use domain parallelization + LAPACK for solving sobspace eigen problem
ifeq ($(USE_DP_SUBEIG), 1)
CPPFLAGS += -DUSE_DP_SUBEIG
endif

# to compile with DEBUG mode
ifeq ($(DEBUG_MODE), 1)
CPPFLAGS += -Wall -g -DDEBUG # -Wno-format-truncation -Wno-stringop-truncation
endif

# to enable SIMD for complex stencil routines
ifeq ($(ENABLE_SIMD_COMPLEX), 1)
CPPFLAGS += -DENABLE_SIMD_COMPLEX
endif

# for old Intel compiler, use -qopenmp instead of -fopenmp. ICC 17 and later also accepts -fopenmp. 
#CFLAGS = -std=gnu99 -O3 -fopenmp
CFLAGS = -std=gnu99 -O3 -fopenmp

MY_OBJS = my_Initialize.o

OBJ_INITIALIZE = \
Initialize.o \
print_usage.o \
SPARC_Input_MPI_create.o \
check_inputs.o \
set_defaults.o \
bcast_SPARC_Atom.o \
SPARC_copy_input.o \
Calculate_SplineDerivRadFun.o \
Cart2nonCart_transformMat.o \
nonCart2Cart_coord.o \
Cart2nonCart_coord.o \
nonCart2Cart_grad.o \
CalculateDistance.o \
write_output_init.o \
Calculate_kpoints.o \
Calculate_local_kpoints.o \
kpointWeight.o \
estimate_memory.o \
Ecut2h.o \
compute_nearest_neighbor_dist.o \
parallel_eigensolver_max_processor.o

OBJLIBS = readfiles.o atomdata.o parallelization.o relax.o tools.o md.o      \
        electrostatics.o electronicGroundState.o electronDensity.o orbitalElecDensInit.o           \
        occupation.o gradVecRoutines.o gradVecRoutinesKpt.o nlocVecRoutines.o     \
        hamiltonianVecRoutines.o lapVecRoutines.o lapVecRoutinesKpt.o \
        linearSolver.o mixing.o exchangeCorrelation.o eigenSolver.o eigenSolverKpt.o energy.o      \
        forces.o stress.o pressure.o finalization.o spinOrbitCoupling.o printing.o \
        linearAlgebra.o \
        xc/vdW/d3/d3correction.o xc/vdW/d3/d3findR0ab.o xc/vdW/d3/d3copyC6.o                       \
        xc/vdW/d3/d3initialization.o xc/vdW/d3/d3finalization.o xc/vdW/d3/d3forceStress.o          \
        xc/vdW/vdWDF/vdWDFinitialization.o xc/vdW/vdWDF/vdWDFfinalization.o                        \
        xc/vdW/vdWDF/vdWDFexchangeLinearCorre.o xc/vdW/vdWDF/vdWDFnonlinearCorre.o                 \
        xc/vdW/vdWDF/vdWDFstress.o xc/vdW/vdWDF/vdWDFparallelization.o                        \
        xc/vdW/vdWDF/vdWDFreadKernel.o \
        xc/mgga/mGGAscan.o xc/mgga/mGGAhamiltonianTerm.o xc/mgga/mGGAstress.o  \
        xc/mgga/mGGAinitialization.o xc/mgga/mGGAfinalization.o xc/mgga/mGGAexchangeCorrelation.o  \
        xc/exx/exactExchange.o xc/exx/exactExchangeKpt.o xc/exx/exactExchangeInitialization.o      \
        xc/exx/exactExchangeFinalization.o xc/exx/exactExchangeStress.o                            \
        xc/exx/exactExchangePressure.o xc/exx/exactExchangeEnergyDensity.o \
        highT/sq.o highT/sqInitialization.o highT/sqFinalization.o highT/sqEnergy.o\
        highT/sqDensity.o highT/sqNlocVecRoutines.o highT/sqParallelization.o          \
        highT/sqProperties.o highT/sqtool.o \
        $(OBJ_INITIALIZE) \
        $(MY_OBJS)

OBJSC = main.o $(OBJLIBS)

EXECUTABLE = sparc.x
TESTBASE = ../.ci

override CC=mpicc

all: lib

# Note the implicit rule to compile '.c' files into '.o' files is
# %.o : %.c
# 	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

lib: $(OBJLIBS)
	rm -fv libsparc.a
	ar rcs libsparc.a $(OBJLIBS)

main_debug: lib
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o main_debug.x main_debug.c libsparc.a $(LDLIBS)

sparc: $(OBJSC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(EXECUTABLE) $^ $(LDLIBS)

.PHONY: clean
clean:
	rm -f  $(OBJSC) $(EXECUTABLE) main_debug.x

test: ../tests/SPARC_testing_script.py
	cd ../tests; python SPARC_testing_script.py
