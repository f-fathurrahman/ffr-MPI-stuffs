static char help[] = "Using DMDA+KSP FD2d to solve Poisson equation";

#include <petsc.h>

PetscErrorCode form_exact_solution(DM da, Vec u_exact)
{
  PetscErrorCode ierr;
  DMDALocalInfo info;
  int i, j;
  double hx, hy, x, y, **au_exact;

  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

  hx = 1.0/(info.mx - 1);
  hy = 1.0/(info.my - 1);

  ierr = DMDAVecGetArray(da, u_exact, &au_exact); CHKERRQ(ierr);
  
  for(j = info.ys; j < info.ys+info.ym; j++) {
    y = j*hy;
    for(i = info.xs; i < info.xs+info.xm; i++) {
      x = i*hx;
      au_exact[j][i] = x*x * (1.0 - x*x) * y*y * (y*y - 1.0);
    }
  }

  ierr = DMDAVecRestoreArray(da, u_exact, &au_exact); CHKERRQ(ierr);

  return ierr;
}


PetscErrorCode 


int main( int argc, char** argv )
{
  PetscErrorCode ierr;

  PetscInitialize( &argc, &argv, (char*)0, help );

  //
  // Initialize DMDA
  //
  DM da;
  ierr = DMDACreate2d( PETSC_COMM_WORLD,
    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
    9,9, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); CHKERRQ(ierr);

  //
  // Linear system matrix
  //
  Mat A;
  //
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMCreateMatrix(da, &A); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);

  //
  // RHS b, approx solution u, exact solution u_exact
  //
  Vec b, u, u_exact;
  //
  ierr = DMCreateGlobalVector(da, &b); CHKERRQ(ierr);
  ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
  ierr = VecDuplicate(b, &u_exact); CHKERRQ(ierr);

  // Fill the vectors
  ierr = form_exact_solution(da, u_exact); CHKERRQ(ierr);
  //ierr = form_RHS(da, b); CHKERRQ(ierr);
  //ierr = form_matrix(da, A); CHKERRQ(ierr); 

  PetscPrintf(PETSC_COMM_WORLD, "Pass here\n");

  PetscFinalize();

  return 0;
}