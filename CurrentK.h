#include "common.h"

#pragma acc routine
real alpha_p(real V);
#pragma acc routine
real beta_p(real V);

#pragma acc routine
real p_inf(real V);

#pragma acc routine
real CurrentK(real V, real p);