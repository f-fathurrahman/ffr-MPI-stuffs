static char HELP_STR[] = "Assemble Mat object using sparsely\n";

#include <petsc.h>

int main( int argc, char** argv )
{
  PetscErrorCode ierr;

  // Row and column indices
  int i1[3] = {0, 1, 2};
  int j1[3] = {0, 1, 2};

  int i2 = 3;
  int j2[3] = {1, 2, 3};

  int i3 = 1;
  int j3 = 3;

  // Values
  double aA1[9] = {1.0, 2.0, 3.0, 
                   2.0, 1.0, -2.0,
                   -1.0, 1.0, 1.0};
  double aA2[3] = {1.0, 1.0, -1.0};
  double aA3 = -3.0;

  ierr = PetscInitialize( &argc, &argv, NULL, NULL ); CHKERRQ( ierr );

  Mat A;
  ierr = MatCreate( PETSC_COMM_WORLD, &A ); CHKERRQ( ierr );
  ierr = MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, 4, 4 ); CHKERRQ( ierr );
  ierr = MatSetFromOptions(A); CHKERRQ( ierr );
  ierr = MatSetUp( A ); CHKERRQ( ierr );

  ierr = MatSetValues( A, 3,i1, 3,j1, aA1, INSERT_VALUES ); CHKERRQ( ierr );
  ierr = MatSetValues( A, 1,&i2, 3,j2, aA2, INSERT_VALUES ); CHKERRQ( ierr );
  ierr = MatSetValue( A, i3, j3, aA3, INSERT_VALUES ); CHKERRQ( ierr );

  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );

  MatDestroy( &A );

  return PetscFinalize();

}
