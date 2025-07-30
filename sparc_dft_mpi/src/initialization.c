/**
 * @file    initialization.c
 * @brief   This file contains the functions for initialization.
 *
 * @authors Qimen Xu <qimenxu@gatech.edu>
 *          Abhiraj Sharma <asharma424@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 *          Hua Huang <huangh223@gatech.edu>
 *          Edmond Chow <echow@cc.gatech.edu>
 *          Alfredo Metere (GPU support), Lawrence Livermore National Laboratory
 * <metere1@llnl.gov>, <alfredo.metere@xsilico.com>
 *
 * Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>
// this is for checking existence of files
#include <unistd.h>
#include "initialization.h"
#include "readfiles.h"
#include "nlocVecRoutines.h"
#include "electrostatics.h"
#include "tools.h"
#include "eigenSolver.h"    // Mesh2ChebDegree, init_GTM_CheFSI()
#include "eigenSolverKpt.h" // init_GTM_CheFSI_kpt()
#include "parallelization.h"
#include "isddft.h"
#include "d3initialization.h"
#include "vdWDFinitialization.h"
#include "mGGAinitialization.h"
#include "exactExchangeInitialization.h"
#include "spinOrbitCoupling.h"
#include "sqInitialization.h"
#include "sqParallelization.h"

#define TEMP_TOL 1e-12

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define N_MEMBR 162

