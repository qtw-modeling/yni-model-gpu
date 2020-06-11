#include "common.h"

#pragma acc routine
real alpha_q(real V);
#pragma acc routine
real beta_q(real V);

#pragma acc routine
real q_inf(real V);

#pragma acc routine
real CurrentHyperpolar(real V, real q);