#include "common.hpp"

// GPU vers of funcs
#pragma acc routine
real alpha_m(real V)
{
    return (V + 37.) / (1. - exp( -(V + 37.)/10. ));
}

// TODO: continue from here and do below!
#pragma acc routine
real beta_m(real V)
{
    return 40. * exp(-5.6 * 1e-2 * (V + 62.));
}

real alpha_h(real V)
{
    return 1.209 * 1e-3 * exp( (-V + 20.)/6.534);
}


real beta_h(real V)
{
    return 1. / ( exp(-(V + 30.)/10.) + 1.);
}


#pragma acc routine
real m_inf(real V)
{
    return alpha_m(V) / (alpha_m(V) + beta_m(V));
}


#pragma acc routine
real h_inf(real V)
{
    return alpha_h(V) / (alpha_h(V) + beta_h(V));
}


#pragma acc routine
real CurrentNa(real V, real m, real n, real h)
{
    //real ENernst = 30.;
    //real gMax = 120.;

    // TODO: check sign;
    // the formula below --- according to the original article of YNI model
    return m * m * m * h * 0.5 * (V - 30.);
}