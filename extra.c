#include "extra.h"


real CalculateLinearInterpolate(real x, real xL, real xR, real yL, real yR)
{
    // calc value at point x \in [xL, xR]
    
    real a = (yR - yL) / (xR - xL); // just the derivative
    real b = yR - a*xR;
    return a*x + b;
}