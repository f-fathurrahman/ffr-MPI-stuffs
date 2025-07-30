## Workflow in Windows

Open WSL2 in terminal, navigate to this directory `src`, then run VSCode `code .`
with `src` as work directory.

Start at `main_debug.c`

Build: `make main_debug`

`my_Initialize`: remove DEBUG (now all debug messages are printed)



## Global variable

It seems that global variable handling all stuffs is the `SPARC_OBJ` struct defined
in `isddft.h`. Most of its members are of primitive types, except for some members
such as:
```c
// ...
PSD_OBJ *psd;           // struct array storing pseudopotential info.
```

It seems to be initialized by calling `Initialize`.

## Input

All input informations seems to be stored 
