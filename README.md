Various Fortran codes utilizing MPI.

OpenMPI version

```
$ mpirun --version
mpirun (Open MPI) 1.10.2
```

Fortran compilers version:

```
ifort (IFORT) 17.0.1 20161005

GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.4) 5.4.0 20160609
```

I could not manage to get `g95` compiler to run with the particular OpenMPI
library that I use.
Error message:
```
In file mpif-sizeof.h:18

    Included at mpif.h:61

    Included at hello.f90:6

        USE, INTRINSIC :: iso_fortran_env, ONLY: REAL128
                                                 1
Error: Symbol 'real128' referenced at (1) not found in module 'iso_fortran_env'
```

