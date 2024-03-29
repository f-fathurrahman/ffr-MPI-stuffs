//STARTWHOLE
static char HELP_STR[] = "Solve a tridiagonal system of arbitrary size.\n"
"Option prefix = tri_.\n";

#include <petsc.h>

int main( int argc, char** argv )
{

  PetscErrorCode ierr;

  PetscInitialize( &argc, &argv, NULL, HELP_STR );

  // Setup problem size, it can be set from program options
  int m = 4;
  ierr = PetscOptionsBegin( PETSC_COMM_WORLD,"tri_", "options for tri", "" ); CHKERRQ( ierr) ;
  ierr = PetscOptionsInt( "-m", "dimension of linear system","tri.c", m, &m, NULL ); CHKERRQ( ierr );
  ierr = PetscOptionsEnd(); CHKERRQ( ierr );

  // Setup vectors
  Vec x, b, x_exact;
    
  ierr = VecCreate( PETSC_COMM_WORLD, &x ); CHKERRQ( ierr );
  ierr = VecSetSizes( x, PETSC_DECIDE, m ); CHKERRQ( ierr );
  ierr = VecSetFromOptions( x ); CHKERRQ( ierr );
    
  ierr = VecDuplicate( x, &b ); CHKERRQ( ierr );
    
  ierr = VecDuplicate( x, &x_exact ); CHKERRQ( ierr );

  // Setup matrix
  Mat A;
  ierr = MatCreate( PETSC_COMM_WORLD, &A ); CHKERRQ( ierr );
  ierr = MatSetSizes( A, PETSC_DECIDE,PETSC_DECIDE, m, m ); CHKERRQ( ierr );
  ierr = MatSetOptionsPrefix( A, "a_"); CHKERRQ( ierr );
  ierr = MatSetFromOptions( A ); CHKERRQ( ierr );
  ierr = MatSetUp( A ); CHKERRQ( ierr );
    
  // Set elements of the matrix and exact (fabricated) solution
  int i, Istart, Iend, j[3];
  double xval;
  double v[3];
  ierr = MatGetOwnershipRange( A, &Istart, &Iend ); CHKERRQ( ierr );
  for( i = Istart; i < Iend; i++ )
  {
    if (i == 0) // first row
    {
      v[0] = 3.0;  v[1] = -1.0;
      j[0] = 0;    j[1] = 1;
      ierr = MatSetValues( A, 1,&i, 2,j, v,INSERT_VALUES); CHKERRQ( ierr );
    }
    else {
      v[0] = -1.0;  v[1] = 3.0;  v[2] = -1.0;
      j[0] = i-1;   j[1] = i;    j[2] = i+1;
      if (i == m-1) // last row: only use first two elements of v and j
      {
        ierr = MatSetValues( A, 1,&i, 2,j, v, INSERT_VALUES ); CHKERRQ( ierr );
      }
      else {
        ierr = MatSetValues( A, 1,&i, 3,j, v, INSERT_VALUES); CHKERRQ( ierr );
      }
    }  
    xval = exp( cos(i) ); // fabricated solution
    ierr = VecSetValues( x_exact, 1,&i, &xval, INSERT_VALUES ); CHKERRQ( ierr );
  }
    
  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
    
  ierr = VecAssemblyBegin( x_exact ); CHKERRQ( ierr );
  ierr = VecAssemblyEnd( x_exact ); CHKERRQ( ierr );
    
  // Set up RHS
  ierr = MatMult( A, x_exact, b ); CHKERRQ( ierr );

  // Setup linear equation solver
  KSP ksp;
  ierr = KSPCreate( PETSC_COMM_WORLD, &ksp ); CHKERRQ( ierr );
  ierr = KSPSetOperators( ksp, A, A ); CHKERRQ( ierr );
  ierr = KSPSetFromOptions( ksp ); CHKERRQ( ierr );
  ierr = KSPSolve( ksp, b, x ); CHKERRQ( ierr );

  // Calculate error
  double errnorm;
  ierr = VecAXPY( x, -1.0, x_exact ); CHKERRQ( ierr );  // x <- x - x_exact
  ierr = VecNorm( x, NORM_2, &errnorm ); CHKERRQ( ierr );
  ierr = PetscPrintf( PETSC_COMM_WORLD,
      "error for m = %d system is |x-x_exact|_2 = %.1e\n",m,errnorm); CHKERRQ( ierr );

  KSPDestroy(&ksp);
  MatDestroy( &A );
  VecDestroy( &x );
  VecDestroy( &b );
  VecDestroy(&x_exact);
    
  return PetscFinalize();
}
