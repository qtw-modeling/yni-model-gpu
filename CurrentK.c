#include "common.h"

// GPU vers of funcs
#pragma acc routine
real alpha_p(real V)
{
    return 9. * 1e-3 / (1. + exp( -(V + 3.8)/9.71)) + 6e-4;
}


#pragma acc routine
real beta_p(real V)
{
    return 2.25 * 1e-4 * (V + 40.) / ( exp((V + 40.) / 13.3) - 1 ) ;
}



#pragma acc routine
real p_inf(real V)
{
    return alpha_p(V) / (alpha_p(V) + beta_p(V));
}


#pragma acc routine
real CurrentK(real V, real p)
{
    // the formulae below --- according to the original article of YNI model
    
    real iK = 0.7 * ( exp(0.0277 * (V + 90.)) - 1.) / exp(0.0277 * (V + 40.));
    
    return p * iK;
}