#include "common.h"

#pragma acc routine
real CurrentLeak(real V)
{
    return 0.8 * ( 1. - exp(-(V + 60.)/20.)  );
}