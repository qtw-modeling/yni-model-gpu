#include "common.hpp"

// GPU vers of funcs
#pragma acc routine
real alpha_q(real V)
{
    return ( 3.4 * 1e-4 * (V + 100.) ) / ( exp((V + 100.)/4.4) - 1.) + 4.95 * 1e-5;
}


#pragma acc routine
real beta_q(real V)
{
    return ( 5 * 1e-4 * (V + 40.) ) / ( 1. - exp(-(V + 40.) / 6.) ) + 8.45 * 1e-5;
}

#pragma acc routine
real q_inf(real V)
{
    return alpha_q(V) / (alpha_q(V) + beta_q(V));
}


real CurrentHyperpolar(real V, real q)
{
    real iH = 0.4 * (V + 25.);
    return q * iH;
}