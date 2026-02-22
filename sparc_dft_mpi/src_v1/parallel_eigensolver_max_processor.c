/**
 * @brief   Simple linear model for selecting maximum number of
 *          processors for eigenvalue solver.
 */
int parallel_eigensolver_max_processor(int N, char RorC, char SorG) {
  if (SorG == 'S') {
    return (int)(0.026918 * N + 1);
  } else if (SorG == 'G') {
    if (RorC == 'R')
      return (int)(0.036215 * N + 1);
    else if (RorC == 'C')
      return (int)(0.038695 * N + 1);
  }
  return -1;
}
