#include "common.h"

#pragma acc routine
real alpha_m(real V);
#pragma acc routine
real beta_m(real V);

#pragma acc routine
real alpha_h(real V);
#pragma acc routine
real beta_h(real V);

#pragma acc routine
real m_inf(real V);
#pragma acc routine
real h_inf(real V);

#pragma acc routine
real CurrentNa(real V, real m, real h);