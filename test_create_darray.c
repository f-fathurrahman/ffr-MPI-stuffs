#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
/* tests */
int darray_2d_c_test1(void);
int darray_4d_c_test1(void);
/* helper functions */
static int parse_args(int argc, char **argv);
static int pack_and_unpack(char *typebuf, int count, MPI_Datatype datatype, int typebufsz);
 
int main(int argc, char **argv)
{
    int err, errs = 0;
 
    MPI_Init(&argc, &argv);
    /* To improve reporting of problems about operations, we change the error handler to errors return */
    MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_RETURN );
 
    /* perform some tests */
    err = darray_2d_c_test1();
    if (err) { fprintf(stderr, "%d errors in 2d darray c test 1.\n", err);fflush(stderr); }
    errs += err;
 
    err = darray_4d_c_test1();
    if (err) { fprintf(stderr, "%d errors in 4d darray c test 1.\n", err);fflush(stderr); }
    errs += err;
 
    MPI_Finalize();
    return errs;
}
 
/* darray_2d_test1()
 *
 * Performs a sequence of tests building darrays with single-element
 * blocks, running through all the various positions that the element might
 * come from.
 *
 * Returns the number of errors encountered.
 */
int darray_2d_c_test1(void)
{
    MPI_Datatype darray;
    int array[9]; /* initialized below */
    int array_size[2] = {3, 3};
    int array_distrib[2] = {MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};
    int array_dargs[2] = {MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG};
    int array_psizes[2] = {3, 3};
    int i, rank, err, errs = 0, sizeoftype;
 
    /* pretend we are each rank, one at a time */
    for (rank=0; rank < 9; rank++) {
        /* set up buffer */
        for (i=0; i < 9; i++) {
            array[i] = i;
        }
        /* set up type */
        err = MPI_Type_create_darray(9 /* size */, rank, 2 /* dims */, array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_INT, &darray);
        if (err != MPI_SUCCESS) {
            errs++;
            fprintf(stderr, "error in MPI_Type_create_darray call; aborting after %d errors\n", errs); fflush(stderr);
            return errs;
        }

        MPI_Type_commit(&darray);
        MPI_Type_size(darray, &sizeoftype);
        if (sizeoftype != sizeof(int)) {
            errs++;
            fprintf(stderr, "size of type = %d; should be %d\n", sizeoftype, (int) sizeof(int));fflush(stderr);
            return errs;
        }

        err = pack_and_unpack((char *) array, 1, darray, 9*sizeof(int));

        for (i=0; i < 9; i++) {
            if ((i == rank) && (array[i] != rank)) {
                errs++;
                fprintf(stderr, "array[%d] = %d; should be %d\n", i, array[i], rank);fflush(stderr);
            }
            else if ((i != rank) && (array[i] != 0)) {
                errs++;
                fprintf(stderr, "array[%d] = %d; should be %d\n", i, array[i], 0);fflush(stderr);
            }
        }
        MPI_Type_free(&darray);
    }
    return errs;
}
 
/* darray_4d_c_test1()
 *
 * Returns the number of errors encountered.
 */
int darray_4d_c_test1(void)
{
    MPI_Datatype darray;
    int array[72];
    int array_size[4] = {6, 3, 2, 2};
    int array_distrib[4] = { MPI_DISTRIBUTE_BLOCK,
                                   MPI_DISTRIBUTE_BLOCK,
                                   MPI_DISTRIBUTE_NONE,
                                   MPI_DISTRIBUTE_NONE };
    int array_dargs[4] = { MPI_DISTRIBUTE_DFLT_DARG,
                                   MPI_DISTRIBUTE_DFLT_DARG,
                                   MPI_DISTRIBUTE_DFLT_DARG,
                                   MPI_DISTRIBUTE_DFLT_DARG };
    int array_psizes[4] = {6, 3, 1, 1};
    int i, rank, err, errs = 0, sizeoftype;
 
    for (rank=0; rank < 18; rank++) {
        /* set up array */
        for (i=0; i < 72; i++) {
            array[i] = i;
        }
 
        /* set up type */
        err = MPI_Type_create_darray(18 /* size */, rank, 4 /* dims */, array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_INT, &darray);
        if (err != MPI_SUCCESS) {
            errs++;
            fprintf(stderr, "error in MPI_Type_create_darray call; aborting after %d errors\n", errs);fflush(stderr);
            return errs;
        }
        MPI_Type_commit(&darray);
        /* verify the size of the type */
        MPI_Type_size(darray, &sizeoftype);
        if (sizeoftype != 4*sizeof(int)) {
            errs++;
            fprintf(stderr, "size of type = %d; should be %d\n", sizeoftype, (int) (4*sizeof(int)));fflush(stderr);
            return errs;
        }
        /* pack and unpack the type, zero'ing out all other values */
        err = pack_and_unpack((char *) array, 1, darray, 72*sizeof(int));
        for (i=0; i < 4*rank; i++) {
            if (array[i] != 0) {
                errs++;
                fprintf(stderr, "array[%d] = %d; should be %d\n", i, array[i], 0);fflush(stderr);
            }
        }
        for (i=4*rank; i < 4*rank + 4; i++) {
            if (array[i] != i) {
                errs++;
                fprintf(stderr, "array[%d] = %d; should be %d\n", i, array[i], i);fflush(stderr);
            }
        }
        for (i=4*rank+4; i < 72; i++) {
            if (array[i] != 0) {
                errs++;
                fprintf(stderr, "array[%d] = %d; should be %d\n", i, array[i], 0);fflush(stderr);
            }
        }
        MPI_Type_free(&darray);
    }
    return errs;
}
 
/******************************************************************/
/* pack_and_unpack()
 *
 * Perform packing and unpacking of a buffer for the purposes of checking
 * to see if we are processing a type correctly. Zeros the buffer between
 * these two operations, so the data described by the type should be in
 * place upon return but all other regions of the buffer should be zero.
 *
 * Parameters:
 * typebuf - pointer to buffer described by datatype and count that
 * will be packed and then unpacked into
 * count, datatype - description of typebuf
 * typebufsz - size of typebuf; used specifically to zero the buffer
 * between the pack and unpack steps
 *
 */
static int pack_and_unpack(char *typebuf, int count, MPI_Datatype datatype, int typebufsz)
{
    char *packbuf;
    int err, errs = 0, pack_size, type_size, position;
 
    err = MPI_Type_size(datatype, &type_size);
    if (err != MPI_SUCCESS) {
        errs++;
        fprintf(stderr, "error in MPI_Type_size call; aborting after %d errors\n", errs);fflush(stderr);
        return errs;
    }
    type_size *= count;
    err = MPI_Pack_size(count, datatype, MPI_COMM_SELF, &pack_size);
    if (err != MPI_SUCCESS) {
        errs++;
        fprintf(stderr, "error in MPI_Pack_size call; aborting after %d errors\n", errs);fflush(stderr);
        return errs;
    }
    packbuf = (char *) malloc(pack_size);
    if (packbuf == NULL) {
        errs++;
        fprintf(stderr, "error in malloc call; aborting after %d errors\n", errs);fflush(stderr);
        return errs;
    }
    position = 0;
    err = MPI_Pack(typebuf, count, datatype, packbuf, type_size, &position, MPI_COMM_SELF);
    if (position != type_size) {
        errs++;
        fprintf(stderr, "position = %d; should be %d (pack)\n", position, type_size);fflush(stderr);
    }
    memset(typebuf, 0, typebufsz);
    position = 0;
    err = MPI_Unpack(packbuf, type_size, &position, typebuf, count, datatype, MPI_COMM_SELF);
    if (err != MPI_SUCCESS) {
        errs++;
        fprintf(stderr, "error in MPI_Unpack call; aborting after %d errors\n", errs);fflush(stderr);
        return errs;
    }
    free(packbuf);
    if (position != type_size) {
        errs++;
        fprintf(stderr, "position = %d; should be %d (unpack)\n", position, type_size);fflush(stderr);
    }
    return errs;
}
