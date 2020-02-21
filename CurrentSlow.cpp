#include "common.hpp"

// GPU vers of funcs
#pragma acc routine
real alpha_d(real V)
{
    return ( 1.045 * 1e-2 * (V + 35.) ) / ( 1. - exp(-(V + 35.)/2.5) ) 
         + ( 3.125 * 1e-2 * V ) / ( 1. - exp(-V/4.8) );
}


#pragma acc routine
real beta_d(real V)
{
    return ( 4.21 * 1e-3 * (V - 5.) ) / ( exp( (V - 5.)/2.5 ) - 1. );
}

real alpha_f(real V)
{
    return ( 3.55 * 1e-4 * (V + 20.) ) / ( exp((V + 20.)/5.633) - 1. );
}


real beta_f(real V)
{
    return ( 9.44 * 1e-4 * fabs(V + 60.) ) / ( 1. + exp(-(V + 29.5)/4.16) );
}


#pragma acc routine
real d_inf(real V)
{
    return alpha_d(V) / (alpha_d(V) + beta_d(V));
}


#pragma acc routine
real f_inf(real V)
{
    return alpha_f(V) / (alpha_f(V) + beta_f(V));
}


#pragma acc routine
real CurrentSlow(real V, real d, real f)
{
    
    // the formula below --- according to the original article of YNI model
    real iS = 12.5 * ( exp((V - 30.)/15.) - 1.); // NB: "(" sign may be misprinted in the YNI-model's article
    
    return ( 0.95*d + 0.05 ) * ( 0.95*f + 0.05) * iS;
}