#include "common.h"

#pragma acc routine
real alpha_d(real V);
#pragma acc routine
real beta_d(real V);

#pragma acc routine
real alpha_f(real V);
#pragma acc routine
real beta_f(real V);

#pragma acc routine
real d_inf(real V);
#pragma acc routine
real f_inf(real V);

#pragma acc routine
real CurrentSlow(real V, real d, real f);