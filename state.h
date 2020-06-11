#pragma once
#include "common.h"

// main (independent) variables 
struct State {
    real V, m, h, p, q, d, f;
};

// extra (dependent) variables
// TODO: add at later modification of the prorgam

struct Currents {
    real* INa, *IK, *ILeak, *IHyperpolar, *ISlow;
};