static char help[] = "Using DMDA+KSP FD2d to solve Poisson equation";

#include <petsc.h>

PetscErrorCode form_exact_solution(DM da, Vec u_exact)
{

  PetscPrintf(PETSC_COMM_WORLD, "\nForming exact solution\n");

  PetscErrorCode ierr;
  DMDALocalInfo info;
  int i, j;
  double hx, hy, x, y, **au_exact;

  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

  hx = 1.0/(info.mx - 1);
  hy = 1.0/(info.my - 1);

  ierr = DMDAVecGetArray(da, u_exact, &au_exact); CHKERRQ(ierr);


  PetscPrintf(PETSC_COMM_WORLD, "info.xs = %d\n", info.xs);
  PetscPrintf(PETSC_COMM_WORLD, "info.xm = %d\n", info.xm);

  PetscPrintf(PETSC_COMM_WORLD, "info.ys = %d\n", info.ys);
  PetscPrintf(PETSC_COMM_WORLD, "info.ym = %d\n", info.ym);

  for(j = info.ys; j < info.ys+info.ym; j++) {
    y = j*hy;
    for(i = info.xs; i < info.xs+info.xm; i++) {
      x = i*hx;
      au_exact[j][i] = x*x * (1.0 - x*x) * y*y * (y*y - 1.0);
    }
  }

  ierr = DMDAVecRestoreArray(da, u_exact, &au_exact); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "Done forming exact solution.\n");

  return ierr;
}


PetscErrorCode form_RHS(DM da, Vec b)
{
  PetscErrorCode ierr;
  int i, j;
  double hx, hy, x, y, f, **ab;
  DMDALocalInfo info;

  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  hx = 1.0/(info.mx-1);
  hy = 1.0/(info.my-1);
  ierr = DMDAVecGetArray(da, b, &ab); CHKERRQ(ierr);
  for( j = info.ys; j < info.ys + info.ym; j++ ) {
    y = j*hy;
    for( i = info.xs; i < info.xs + info.xm; i++ ) {
      x = i*hx;
      if( i==0 || i==info.mx-1 || j==0 || j==info.my-1 ) {
        ab[j][i] == 0.0; // on boundary: 1*u = 0
      }
      else {
        f = 2.0*( (1.0 - 6.0*x*x) * y*y*(1.0 - y*y) +
                  (1.0 - 6.0*y*y) * x*x*(1.0 - x*x) );
        ab[j][i] = hx * hy * f;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da, b, &ab); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode form_matrix( DM da, Mat A )
{

  PetscPrintf(PETSC_COMM_WORLD, "\nForming matrix.\n");

  PetscErrorCode ierr;
  DMDALocalInfo info;
  MatStencil row, col[5];
  double hx, hy, v[5];
  int i, j, Ncols;

  int ic;

  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  //
  hx = 1.0/(info.mx - 1);
  hy = 1.0/(info.my - 1);
  //
  PetscPrintf(PETSC_COMM_WORLD, "hx = %g\n", hx);
  PetscPrintf(PETSC_COMM_WORLD, "hy = %g\n", hy);
  //
  PetscPrintf(PETSC_COMM_WORLD, "hx/hy = %g\n", hx/hy);
  PetscPrintf(PETSC_COMM_WORLD, "hy/hx = %g\n", hy/hx);
  //
  double hx_hy = hx/hy;
  double hy_hx = hy/hx;

  //
  for( j = info.ys; j < info.ys+info.ym; j++ ) {
    for( i = info.xs; i < info.xs+info.xm; i++ ) {
      //
      row.j = j;  // row of A corresponding to (x_i,y_j)
      row.i = i;
      //
      col[0].j = j;
      col[0].i = i;
      Ncols = 1;
      //
      if( i==0 || i==info.mx-1 || j==0 || j==info.my-1 ) {
        v[0] = 1.0; // on boundary
      }
      else {
        //
        PetscPrintf(PETSC_COMM_WORLD, "Pass here 0\n");
        v[0] = 2*(hy_hx + hx_hy); // interior
        //
        if( i-1 > 0 ) {
          PetscPrintf(PETSC_COMM_WORLD, "Pass here 1\n");
          col[Ncols].j = j;
          col[Ncols].i = i-1;
          v[Ncols++] = -hy_hx;
        }
        if( i+1 < info.mx-1 ) {
          PetscPrintf(PETSC_COMM_WORLD, "Pass here 2\n");
          col[Ncols].j = j;
          col[Ncols].i = i+1;
          v[Ncols++] = -hy_hx;
        }
        if( j-1 > 0 ) {
          PetscPrintf(PETSC_COMM_WORLD, "Pass here 3\n");
          col[Ncols].j = j-1;
          col[Ncols].i = i;
          v[Ncols++] = -hx_hy;
        }
        if( j+1 < info.my-1 ) {
          PetscPrintf(PETSC_COMM_WORLD, "Pass here 4\n");
          col[Ncols].j = j+1;
          col[Ncols].i = i;
          v[Ncols++] = -hx_hy;
        }
      }
      PetscPrintf(PETSC_COMM_WORLD, "Ncols = %d\n", Ncols);
      PetscPrintf(PETSC_COMM_WORLD, "Setting stencil values (global index): %d %d\n", i, j);
      PetscPrintf(PETSC_COMM_WORLD, "Column index:\n");
      for(ic = 0; ic < 5; ic++) {
        PetscPrintf(PETSC_COMM_WORLD, "v = %18.10f ", v[ic]);
        PetscPrintf(PETSC_COMM_WORLD, "[%d %d %d]\n", ic, col[ic].i, col[ic].j);
      }
      ierr = MatSetValuesStencil(A, 1, &row, Ncols, col, v, INSERT_VALUES); CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "\nDone forming matrix.\n");

  return 0;
}


int main( int argc, char** argv )
{
  PetscErrorCode ierr;

  PetscInitialize(&argc, &argv, (char*)0, help);

  PetscMPIInt rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscPrintf(PETSC_COMM_WORLD, "Number of processors = %d, rank = %d\n", size, rank);

  //
  // Initialize DMDA
  //
  // Default size: 9x9
  // Change using option: -da_grid_x M -da_grid_y N
  DM da;
  ierr = DMDACreate2d( PETSC_COMM_WORLD,
    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
    9, 9, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); CHKERRQ(ierr);

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
  ierr = form_RHS(da, b); CHKERRQ(ierr);
  ierr = form_matrix(da, A); CHKERRQ(ierr);

  // Create and solve the linear system
  KSP ksp;
  //
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, u); CHKERRQ(ierr);

  // Report on grid and numerical error
  ierr = VecAXPY(u, -1.0, u_exact); CHKERRQ(ierr); // u <- u + (-1.0)*u_exact
  //
  double errnorm;
  ierr = VecNorm(u, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  //
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  //
  ierr = PetscPrintf(PETSC_COMM_WORLD, "On %d x and %d grid: error |u-u_exact|_inf = %g\n",
    info.mx, info.my, errnorm); CHKERRQ(ierr);


  VecDestroy(&u);
  VecDestroy(&u_exact);
  VecDestroy(&b);

  MatDestroy(&A);
  KSPDestroy(&ksp);
  DMDestroy(&da);

  PetscFinalize();

  return 0;
}