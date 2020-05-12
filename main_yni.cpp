//
// Created by QuoTheWhite on 27/03/2019.
//
//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <fstream>
//#include <cstdlib>
//#include "openacc.h"


#include <string>
#include <map>
#include <cmath>
#include <math.h> // 4some specific constants

// 4 C++11 includes
//#include <vector>

//typedef double real;
#include "state.hpp"
#include "common.hpp" // one for all includes
#include "CurrentNa.hpp"
#include "CurrentK.hpp"
#include "CurrentLeak.hpp"
#include "CurrentHyperpolar.hpp"
#include "CurrentSlow.hpp"

// grid parameters
#define numSegmentsX 30 //10
#define numSegmentsY 30 //10
#define numPointsX (numSegmentsX + 1) // = numCells
#define numPointsY (numSegmentsY + 1) // = numCells
#define numPointsTotal (numPointsX * numPointsY)
#define hx 0.07 //1. // uncomment if cells are connected // (1./numSegmentsX)
#define hy 0.07 //1. // uncomment if cells are connected // (1./numSegmentsY)
#define T (7000) //(1000.) // old val: 500 // endtime (in ms)
#define dt 0.005 // old val = 1e-4 // timestep (in ms)

// model parameters
#define Cm 1.
#define VRest (-60.) // NOTE: there exists no resting potential for SA node

// tissue parameters
#define Dx (7e-3) //88e-3 //1e-3 //1e-3 // conductivity
#define Dy (7e-3) //88e-3 //1e-3 //1e-3 // conductivity



void Write2VTK(std::string fileName, const int n, real* p, const real h, const int step)
{
    // C style
    //char fn[256];
    //sprintf(fn, "./output/yni.%d.vtk", step);

    // C++ style
    std::string fn = "./output/" + fileName;
    char fnEnding[256];
    sprintf(fnEnding, ".%d.vtk", step);
    std::string str_fnEnding = fnEnding;
    fn += str_fnEnding;

    std::fstream f(fn, std::ios::out);
    f << "# vtk DataFile Version 3.0" << std::endl;
    f << "Solution" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET RECTILINEAR_GRID" << std::endl;
    f << "DIMENSIONS " << n + 1 << " " << n + 1 << " 1" << std::endl;
    f << "X_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
    //for (int i = 1; i < n; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Y_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
    //for (int i = 1; i < n; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Z_COORDINATES 1 double\n0" << std::endl;
    f << "CELL_DATA " << (n * n) << std::endl;
    //f << "CELL_DATA " << (n-2) * (n-2) << std::endl;
    f << "SCALARS V_membrane double\nLOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++)
            f << p[j * n + i] << " ";
        f << std::endl;
    }
    f.close();
}


int CalculateLinearCoordinate_CPU(int i, int j) {
    return i + j*numPointsX;
}

#pragma acc routine
int CalculateLinearCoordinate(int i, int j) {
    return i + j*numPointsX;
}


// gating vars kinetic functions' definitions
//////////////////////////////////////////////

// CPU versions of functions

real alpha_n_CPU(real V) {
    return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
}


real beta_n_CPU(real V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}


real n_inf_CPU(real V) {
    return alpha_n_CPU(V) / (alpha_n_CPU(V) + beta_n_CPU(V));
}


real alpha_m_CPU(real V) {
    return 0.1*(V + 40.0)/(1.0 - exp(-(V + 40.0)/10.0));
}


real beta_m_CPU(real V) {
    return 4.0*exp(-(V + 65.0) / 18.0);
}


real m_inf_CPU(real V) {
    return alpha_m_CPU(V) / (alpha_m_CPU(V) + beta_m_CPU(V));
}


real alpha_h_CPU(real V) {
    return 0.07*exp(-(V + 65.5) / 20.0);
}


real beta_h_CPU(real V) {
    return 1.0 / (1.0 + exp(-(V + 35.5) / 10.0));
}


real h_inf_CPU(real V) {
    return alpha_h_CPU(V) / (alpha_h_CPU(V) + beta_h_CPU(V));
}


// GPU versions of functions
/*
#pragma acc routine
real alpha_n(real V) {
    return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
}

#pragma acc routine
real beta_n(real V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}

#pragma acc routine
real n_inf(real V) {
    return alpha_n(V) / (alpha_n(V) + beta_n(V));
}
*/
/*
#pragma acc routine
real alpha_m(real V) {
    return 0.1*(V + 40.0)/(1.0 - exp(-(V + 40.0)/10.0));
}

#pragma acc routine
real beta_m(real V) {
    return 4.0*exp(-(V + 65.0) / 18.0);
}

#pragma acc routine
real m_inf(real V) {
    return alpha_m(V) / (alpha_m(V) + beta_m(V));
}

#pragma acc routine
real alpha_h(real V) {
    return 0.07*exp(-(V + 65.5) / 20.0);
}

#pragma acc routine
real beta_h(real V) {
    return 1.0 / (1.0 + exp(-(V + 35.5) / 10.0));
}

#pragma acc routine
real h_inf(real V) {
    return alpha_h(V) / (alpha_h(V) + beta_h(V));
}
*/

#pragma acc routine
real I_Stim(int i, int j, real value)
{
    //int x0 = (int)((real)numSegmentsX /2.); int y0 = (int)((real)numSegmentsY /2.); // tmp values

    /* uncomment if cells are connected
    if ((i > 1 && i < numPointsX - 2) && (j == 2))
        return value;
    else
        return 0.;
    */

    return value;
}

#pragma acc routine
real TotalIonCurrent(int idx, real V, real m, /* real n,*/ real h, real p, real q, real d, real f,
                     real *INa, real *IK, real *ILeak, real *IHyperpolar, real *ISlow)
{
    INa[idx] = CurrentNa(V, m, h);
    IK[idx] = CurrentK(V, p);
    ILeak[idx] = CurrentLeak(V);
    IHyperpolar[idx] = CurrentHyperpolar(V, q);
    ISlow[idx] = CurrentSlow(V, d, f);

    // TODO: check the sign of the expression below
    return -(INa[idx] + IK[idx] + ILeak[idx] + IHyperpolar[idx] + ISlow[idx]);
}

real TimeViaPhase(real phase_, real Period_, real t0_) 
{
    /* phi = phase; Period = T; t0 = time of first depolarization. */
    
    return t0_ + Period_ / (2.*M_PI) * phase_;
    
}

// phase calculations
/*real*/ std::map<std::string, real> VviaPhase(real phase) 
{
    // we dont need to perform memalloc for members within these structs
    State* Old = new State;
    State* New = new State;
    // initial conditions: only for Old structs: New will be calculated in the loop
    Old->V = VRest;
    Old->m = 0.067;                //m_inf_CPU(VRest); // 0.5
    Old->h = 0.999; //h_inf_CPU(VRest); // 0.5
    Old->p = 0.;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->q = 0.;    // may be false; TODO perform calcs with higher T
    Old->d = 0.; //0.;
    Old->f = 1.; //1.;


    // ...but here --- YES
    Currents* CurrentsOld = new Currents;
    
    // ??? maybe we dont need CurrentsNew??? we are not swapping OLD and NEW
    Currents* CurrentsNew = new Currents;
    // dont forget to memalloc pointers within Currents struct;
    // we use only single cell
    CurrentsOld->INa = new real[1];
    CurrentsOld->IK = new real[1];
    CurrentsOld->ILeak = new real[1];
    CurrentsOld->IHyperpolar = new real[1];
    CurrentsOld->ISlow = new real[1];
    // ...same for "New" struct here
    CurrentsNew->INa = new real[1];
    CurrentsNew->IK = new real[1];
    CurrentsNew->ILeak = new real[1];
    CurrentsNew->IHyperpolar = new real[1];
    CurrentsNew->ISlow = new real[1];

    real THRESHOLD = -30; // hardcoded for now
    real time0; // random value
    real time1; // random value
    bool isThresholdFound = false;

    real tCurrent = 0;
    int counter = 0;
    // TODO:
    // main loop: timestepping
    while (1)
    {

        ///////////////// gating variables: ode ("reaction") step
        // TODO: make only ONE read of Old->V, etc. from memory; to more speedup, esp. for GPU
        New->m = m_inf(Old->V) + (Old->m - m_inf(Old->V)) * exp(-dt * (alpha_m(Old->V) + beta_m(Old->V)));

                    
        New->h = h_inf(Old->V) + (Old->h - h_inf(Old->V)) * exp(-dt * (alpha_h(Old->V) + beta_h(Old->V)));

        New->p = p_inf(Old->V) + (Old->p - p_inf(Old->V)) * exp(-dt * (alpha_p(Old->V) + beta_p(Old->V)));

        New->q = q_inf(Old->V) + (Old->q - q_inf(Old->V)) * exp(-dt * (alpha_q(Old->V) + beta_q(Old->V)));

        New->d = d_inf(Old->V) + (Old->d - d_inf(Old->V)) * exp(-dt * (alpha_d(Old->V) + beta_d(Old->V)));

        New->f = f_inf(Old->V) + (Old->f - f_inf(Old->V)) * exp(-dt * (alpha_f(Old->V) + beta_f(Old->V)));
                  
        // membrane potential calc                  
        New->V = Old->V + dt / Cm * ( TotalIonCurrent(0 /* 0 --- index, stands for the single cell model */ , 
                                 Old->V, Old->m, Old->h, Old->p, Old->q, Old->d, Old->f,
                                 (CurrentsOld->INa), (CurrentsOld->IK), (CurrentsOld->ILeak), 
                                 (CurrentsOld->IHyperpolar), 
                                 (CurrentsOld->ISlow))
                                 + I_Stim(0, 0, 0.*1e0) ); // "standart" I_stim = 1e

        // when threshold is found
        if ((New->V > THRESHOLD) && (Old->V < THRESHOLD))
        {
            // when 2nd threshold time (t1) is found: set t1 and then exit the loop
            if (isThresholdFound == true)
            {

                time1 = tCurrent; // nearest-neighbour interpolaion; change to linear!
                break;            // phase(V)
            }
            else // when threshold time (t0) is found: set t0
            {
                time0 = tCurrent; // nearest-neighbour interpolaion; change to linear!
                isThresholdFound = true;
                //return ; // phase(V)
            }

        }

        
        tCurrent += dt;
        counter += 1;

        // swapping time layers
        State* tmp;
        tmp = Old; Old = New; New = tmp;

        //printf("Iteration #%d\n", counter);

    } // while

    //printf("t0 = %.2f, t1 = %.2f\n", time0, time1);
    //std::cin.get();

    // set vars, calculated within the loop
    real period = (time1 - time0);//*0.5; // period of oscillations; remove "0.5" when period calc bug is found!
    //printf("First loop is finished; period of oscillations: %.2f ms\n", period);
    
    
    // repeat the loop (calculations) again and find V(phi)
    tCurrent = 0; // again
    
    real tOfPhase = TimeViaPhase(phase, period, time0);
    //printf("Phase: %.2f, tOfPhase: %.2f\n", phase, tOfPhase);
    //std::cin.get();

    //real VOfPhase; // to be determined in the loop below
    std::map<std::string, real> stateOfPhase;

    // (again) initial conditions: only for Old structs: New will be calculated in the loop
    Old->V = VRest;
    Old->m = 0.067; //m_inf_CPU(VRest); // 0.5
    Old->h = 0.999; //h_inf_CPU(VRest); // 0.5
    Old->p = 0.;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->q = 0.;    // may be false; TODO perform calcs with higher T
    Old->d = 0.;    //0.;
    Old->f = 1.;    //1.;

    // (again): main loop: timestepping
    while (1)
    {
        // it means, dat we found the moment of time, corresponding to the phase value
        if (tCurrent >= tOfPhase)
        {    
            //VOfPhase = Old->V; // nearest-neighbour iterpolation; change to linear!
            stateOfPhase["V"] = Old->V;
            stateOfPhase["m"] = Old->m;
            stateOfPhase["h"] = Old->h;
            stateOfPhase["p"] = Old->p;
            stateOfPhase["q"] = Old->q;
            stateOfPhase["d"] = Old->d;
            stateOfPhase["f"] = Old->f;
            //stateOfPhase[""]
            //stateOfPhase[""]


                break;
            //return VOfPhase;
        }

        ///////////////// gating variables: ode ("reaction") step
        // TODO: make only ONE read of Old->V, etc. from memory; to more speedup, esp. for GPU
        New->m = m_inf(Old->V) + (Old->m - m_inf(Old->V)) * exp(-dt * (alpha_m(Old->V) + beta_m(Old->V)));

        New->h = h_inf(Old->V) + (Old->h - h_inf(Old->V)) * exp(-dt * (alpha_h(Old->V) + beta_h(Old->V)));

        New->p = p_inf(Old->V) + (Old->p - p_inf(Old->V)) * exp(-dt * (alpha_p(Old->V) + beta_p(Old->V)));

        New->q = q_inf(Old->V) + (Old->q - q_inf(Old->V)) * exp(-dt * (alpha_q(Old->V) + beta_q(Old->V)));

        New->d = d_inf(Old->V) + (Old->d - d_inf(Old->V)) * exp(-dt * (alpha_d(Old->V) + beta_d(Old->V)));

        New->f = f_inf(Old->V) + (Old->f - f_inf(Old->V)) * exp(-dt * (alpha_f(Old->V) + beta_f(Old->V)));

        // membrane potential calc
        New->V = Old->V + dt / Cm * ( TotalIonCurrent(0 /* 0 --- index, stands for the single cell model */ , 
                                 Old->V, Old->m, Old->h, Old->p, Old->q, Old->d, Old->f,
                                 (CurrentsOld->INa), (CurrentsOld->IK), (CurrentsOld->ILeak), 
                                 (CurrentsOld->IHyperpolar), 
                                 (CurrentsOld->ISlow)) 
                                 + I_Stim(0, 0, 0. * 1e0) ); // "standart" I_stim = 1e


        tCurrent += dt;
        //stepNumber += 1;

        // swapping time layers
        State* tmp;
        tmp = Old; Old = New; New = tmp;

    } // while

    //printf("Second loop is finished; VOfPhase: %.1f mV\n", VOfPhase);
    //std::cin.get();
    // "return" --- is within the loop (look up)
    return stateOfPhase; //VOfPhase;
}


void SetInitialConditions_CPU(real* V, real* m, /* real* n,*/ real* h, real* p, real* q, real* d, 
real* f, real value) {
    int idx;
    std::srand(unsigned(1.)); // initial seed for random number generator
    real randomNumber;

    // single initial peak
    //int iCenter = (int)((real)numSegmentsX /2.);
    //int jCenter = (int)((real)numSegmentsY /2.); // tmp values
    //int idxCenter = CalculateLinearCoordinate_CPU(iCenter, jCenter);

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) {

            int idxCenter = CalculateLinearCoordinate_CPU(i, j);
            randomNumber =  ((real)(std::rand() % 20))/20.; // 4phase setting

            // the borders: Dirichlet boundary conditions
            //if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1)) {
                // TODO: find out about the values
                
                // TODO: wrong formulas below ////////////////////////////////////////////////////
                //int jTilde = numPointsY/2 - j; // j = j^{star}; jTilde = jTilde_{star}
                //int iTilde = i - (numPointsX/2); // i = i^{star}; iTilde = iTilde_{star}

                // for phase calculation: using angle in polar coords
                real LTotal = numPointsX*hx; // should = numPointsY*hy
                real L = (j*hy + hy/2.) - (numPointsY*hy/2.); // LTotal/2. - (j*hx + hx/2.) --- old formula, when j-order was incorrect
                real lsmall = LTotal/2. - ( (numPointsX - 1 - i)*hx + hx/2.);

                real phase = atan2(L, lsmall); // = angle in polar coords; use atan2() func instead of atan() !
                
                // check sign: atan2() returns vals. from [-pi, pi]
                if (phase < 0)
                {
                    phase += 2*M_PI;
                }
                

                //real phaseShifted = phase - M_PI/2.; // phase from R.Syunyaev article
                
                //printf("Phase = %.2f deg.\n", phase*180/M_PI);
                //std::cin.get();
                // TODO //////////////////////////////////////////////////////////////

                // the func returns a std::map of all the vars' values
                std::map<std::string, real> stateForPhase = VviaPhase(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase["V"]);
                //std::cin.get();

                V[idxCenter] = stateForPhase["V"];  //VviaPhase(phase); //M_PI/12. //VRest;
                m[idxCenter] = stateForPhase["m"]; //0.067;//m_inf_CPU(VRest); // 0.5
                h[idxCenter] = stateForPhase["h"]; //0.999; //h_inf_CPU(VRest); // 0.5
                p[idxCenter] = stateForPhase["p"]; //0.;// 0.1; // may be false; TODO perform calcs with higher T 
                q[idxCenter] = stateForPhase["q"]; //0.; // may be false; TODO perform calcs with higher T  

                d[idxCenter] = stateForPhase["d"]; //0.; //0.;
                f[idxCenter] = stateForPhase["f"]; //1.; //1.;

                // for progress checking: in percents
                printf("Set. initial cond: %.2f percent completed\n", 
                        100.*idxCenter / CalculateLinearCoordinate_CPU(numSegmentsX, numSegmentsY));
            }

    // after filling the whole area: "fill" borders wiht Neumann boundary cond.
    // the borders: Neumann boundary conditions
    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++)
        {
            int idxCenter = CalculateLinearCoordinate_CPU(i, j);
            
            // borrder cells, including corner cells
            if (i == 0 || j == 0 || i == (numSegmentsX) || j == (numSegmentsY))
            {
                int idxNear;

                if ((i == 0)) //&& (j >= 1) && (j <= numSegmentsY - 1)) // left border, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i + 1, j);
                if ((j == 0)) //&& (i >= 1) && (i <= numSegmentsX - 1)) // bottom, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j + 1);
                if ((j == numSegmentsY)) // && (i >= 1) && (i <= numSegmentsX - 1)) // top, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j - 1);
                if ((i == numSegmentsX)) // && (j >= 1) && (j <= numSegmentsY - 1)) // right, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i - 1, j);

                // what about corner cells? for now, they are not treated (?)
                V[idxCenter] = V[idxNear];
                m[idxCenter] = m[idxNear];
                h[idxCenter] = h[idxNear];
                p[idxCenter] = p[idxNear];
                q[idxCenter] = q[idxNear];
                d[idxCenter] = d[idxNear];
                f[idxCenter] = f[idxNear];
            }
        }
}


/*
#pragma acc routine
real I_Stim(int i, int j, real value) {
    //int x0 = (int)((real)numSegmentsX /2.); int y0 = (int)((real)numSegmentsY /2.); // tmp values

    // // uncomment if cells are connected
    //if ((i > 1 && i < numPointsX - 2) && (j == 2))
    //  return value;
    //else
    //    return 0.;
    ///

   return value;
} */


// ion currents
///////////////////////////////////////////////////
/*
#pragma acc routine
real CurrentNa(real V, real m, real n, real h) {
    real ENernst = 50.;
    real gMax = 120.;

    // TODO: check sign
    return gMax*m*m*m*h*(V - ENernst);
}

#pragma acc routine
real CurrentK(real V, real m, real n, real h) {
    real ENernst = -77.;
    real gMax = 36.;

    // TODO: check sign
    return gMax*n*n*n*n*(V - ENernst);
}

#pragma acc routine
real CurrentLeak(real V, real m, real n, real h) {
    real ENernst = -54.387;
    real gMax = 0.3;

    return gMax*(V - ENernst);
}
*/

/*
#pragma acc routine
real TotalIonCurrent(int idx, real V, real m, real h, real p, real q, real d, real f,
real* INa, real* IK, real* ILeak, real* IHyperpolar, real* ISlow) {
    INa[idx] = CurrentNa(V, m, h);
    IK[idx] = CurrentK(V, p);
    ILeak[idx] = CurrentLeak(V);
    IHyperpolar[idx] = CurrentHyperpolar(V, q);
    ISlow[idx] = CurrentSlow(V, d, f);

    // TODO: check the sign of the expression below
    return  -(INa[idx] + IK[idx] + ILeak[idx] + IHyperpolar[idx] + ISlow[idx]);
}
*/



int main() {

    // setting a GPU for the computations; NOTE: req "openacc.h"!
    //acc_set_device_num(1, acc_device_nvidia);

    // allocating memory
    
    // C++11 style alloc
    //std::vector<real> vecVOld(numPointsTotal);
    //real* VOld = &vecVOld.front();

    // gating vars
    real* VOld = new real[numPointsTotal];
    real* mOld = new real[numPointsTotal];
    //real* nOld = new real[numPointsTotal];
    real* hOld = new real[numPointsTotal];
    real* pOld = new real[numPointsTotal];
    real* qOld = new real[numPointsTotal];
    real* dOld = new real[numPointsTotal];
    real* fOld = new real[numPointsTotal];

    real* VNew = new real[numPointsTotal];
    real* mNew = new real[numPointsTotal];
    //real* nNew = new real[numPointsTotal];
    real* hNew = new real[numPointsTotal];
    real* pNew = new real[numPointsTotal];
    real* qNew = new real[numPointsTotal];
    real* dNew = new real[numPointsTotal];
    real* fNew = new real[numPointsTotal];

    // currents
    real* INaOld = new real[numPointsTotal];
    real* IKOld = new real[numPointsTotal];
    real* ILeakOld = new real[numPointsTotal];
    real* IHyperpolarOld = new real[numPointsTotal];
    real* ISlowOld = new real[numPointsTotal];

    real* INaNew = new real[numPointsTotal];
    real* IKNew = new real[numPointsTotal];
    real* ILeakNew = new real[numPointsTotal];
    real* IHyperpolarNew = new real[numPointsTotal];
    real* ISlowNew = new real[numPointsTotal];

    real* tmp; // a pointer for swapping time-layers 'n' and 'n+1'

    // for output in a loop
    std::map<std::string, real*> variables;
    variables["V"] = VOld;
    variables["m"] = mOld;
    variables["h"] = hOld;
    variables["p"] = pOld;
    variables["q"] = qOld;
    variables["d"] = dOld;
    variables["f"] = fOld;
    variables["INa"] = INaOld;
    variables["IK"] = IKOld;
    variables["ILeak"] = ILeakOld;
    variables["IHyperpolar"] = IHyperpolarOld;
    variables["ISlow"] = ISlowOld;
    // = {"V", "m", "h", "p", "q", "d", "f"};


    // initializing before timesteppin'
    SetInitialConditions_CPU(VOld, mOld, /* nOld,*/ hOld, pOld, qOld, dOld, fOld, 0.);
    //SetInitialConditions_CPU(VNew, mNew, /* nNew,*/ hNew, pNew, qNew, dNew, fNew, 0.); // for avoiding "junk" values in all '...New' arrays

    real tCurrent = 0.;
    int stepNumber = 0;
    int counterOutput = 1;

    printf("Timesteppin' begins...\n");
    clock_t start = clock();

// pragmas without "-acc" flag --- are ignored?
#pragma acc data copy(VOld[0:numPointsTotal], mOld[0:numPointsTotal], nOld[0:numPointsTotal], hOld[0:numPointsTotal], \
		      VNew[0:numPointsTotal], mNew[0:numPointsTotal], nNew[0:numPointsTotal], hNew[0:numPointsTotal]) \
		      deviceptr(tmp)
{
    // main loop: timestepping
    while (tCurrent < T) 
    {

        // TODO: change order of indexing (i, j)
        
	#pragma acc kernels \
	present(VOld[0:numPointsTotal], mOld[0:numPointsTotal], nOld[0:numPointsTotal], hOld[0:numPointsTotal], \
                VNew[0:numPointsTotal], mNew[0:numPointsTotal], nNew[0:numPointsTotal], hNew[0:numPointsTotal])
	{
	
	#pragma acc loop collapse(2) independent
	for (int j = 0; j < numPointsY; j++)
            for (int i = 0; i < numPointsX; i++) 
            {

                int idxCenter = CalculateLinearCoordinate(i, j);
                
                // inner cells
                if (i >= 1 && j >= 1 && i <= (numSegmentsX - 1) && j <= (numSegmentsY - 1))
                {
                    // for short names
                    int idxUp = CalculateLinearCoordinate(i, j + 1);
                    int idxDown = CalculateLinearCoordinate(i, j - 1);
                    int idxLeft = CalculateLinearCoordinate(i - 1, j);
                    int idxRight = CalculateLinearCoordinate(i + 1, j);

                    
                    ///////////////// gating variables: ode ("reaction") step

                    // TODO: make only ONE read of VOld[idxCenter], etc from memory; to more speedup, esp. for GPU
                    mNew[idxCenter] = m_inf(VOld[idxCenter]) + (mOld[idxCenter] - m_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_m(VOld[idxCenter]) + beta_m(VOld[idxCenter])));

                    //nNew[idxCenter] = n_inf(VOld[idxCenter]) + (nOld[idxCenter] - n_inf(VOld[idxCenter]))
                    //                                            * exp(-dt * (alpha_n(VOld[idxCenter]) + beta_n(VOld[idxCenter])));

                    hNew[idxCenter] = h_inf(VOld[idxCenter]) + (hOld[idxCenter] - h_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_h(VOld[idxCenter]) + beta_h(VOld[idxCenter])));
                    
                    pNew[idxCenter] = p_inf(VOld[idxCenter]) + (pOld[idxCenter] - p_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_p(VOld[idxCenter]) + beta_p(VOld[idxCenter])));
                   
                    qNew[idxCenter] = q_inf(VOld[idxCenter]) + (qOld[idxCenter] - q_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_q(VOld[idxCenter]) + beta_q(VOld[idxCenter])));
                    

                    dNew[idxCenter] = d_inf(VOld[idxCenter]) + (dOld[idxCenter] - d_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_d(VOld[idxCenter]) + beta_d(VOld[idxCenter])));
                    
                    fNew[idxCenter] = f_inf(VOld[idxCenter]) + (fOld[idxCenter] - f_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_f(VOld[idxCenter]) + beta_f(VOld[idxCenter])));
                   
                    /*
                    // for steady state's calculation
                    VNew[idxCenter] = VRest;
                    */

                    
                    //////////////////
                    // "discrete diffusion" step
                    VNew[idxCenter] = VOld[idxCenter]
                    // uncomment if cells are connected; otherwize --- Nans
                     + dt / Cm * (
                            Dx  * (VOld[idxRight] - 2 * VOld[idxCenter] + VOld[idxLeft])
                            + Dy  * (VOld[idxUp] - 2 * VOld[idxCenter] + VOld[idxDown])
                    );
                    
                    
                    // reaction step
                    VNew[idxCenter] += dt / Cm * (TotalIonCurrent(idxCenter, VOld[idxCenter], mOld[idxCenter],
                                                             hOld[idxCenter], pOld[idxCenter], 
                                                            qOld[idxCenter], dOld[idxCenter], fOld[idxCenter],
                                                            INaOld, IKOld, ILeakOld, IHyperpolarOld, ISlowOld)
                                                                        + I_Stim(i, j, 0.*1e0)); // "standart" I_stim = 1e0;
                    
               } // if
               
               // the borders: Neumann boundary conditions
               else
               {
                    int idxNear;
                    
                    if ((i == 0) && (j >= 1) && (j <= numSegmentsY - 1)) // left border, except for corner cells
                        idxNear = CalculateLinearCoordinate(i + 1, j);
                    else if ((j == 0) && (i >= 1) && (i <= numSegmentsX - 1)) // bottom, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j + 1);
                    else if ((j == numSegmentsY) && (i >= 1) && (i <= numSegmentsX - 1)) // top, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j - 1);
                    else if ((i == numSegmentsX) && (j >= 1) && (j <= numSegmentsY - 1)) // right, except for corner cells
                        idxNear = CalculateLinearCoordinate(i - 1, j);
                    else { // if corner cell
                        continue; // do nothing, continue the "i,j" loop
                    }

                    // what about corner cells? for now, they are not treated (?)
                    // Neumann boundary cond setting
                    VNew[idxCenter] = VNew[idxNear];
                    mNew[idxCenter] = mNew[idxNear];
                    hNew[idxCenter] = hNew[idxNear];
                    pNew[idxCenter] = pNew[idxNear];
                    qNew[idxCenter] = qNew[idxNear];
                    dNew[idxCenter] = dNew[idxNear];
                    fNew[idxCenter] = fNew[idxNear];
               }

            } // for
	
	} // acc kernels

        if ( (stepNumber % 2000) == 0) { // output each 10 msec: 10/dt = 2000
        //if ( (stepNumber) % (int)(T/dt/500)  == 0 ) {
        #pragma acc update host(VOld[0:numPointsTotal])
            variables["V"] = VOld;
            variables["m"] = mOld;
            variables["h"] = hOld;
            variables["p"] = pOld;
            variables["q"] = qOld;
            variables["d"] = dOld;
            variables["f"] = fOld;
            variables["INa"] = INaOld;
            variables["IK"] = IKOld;
            variables["ILeak"] = ILeakOld;
            variables["IHyperpolar"] = IHyperpolarOld;
            variables["ISlow"] = ISlowOld;

            for(const auto& variable: variables) {// variables repr "X"Old values
                int outNumber = stepNumber;
                
                if ((variable.first.compare("V")) == 0) // output only "V"
                {
                Write2VTK(variable.first, numPointsX, variable.second, hx, outNumber); // for now: numPointsX == numPointsY
                //Write2VTK("V", numPointsX, variables["V"], hx, counterOutput); // for now: numPointsX == numPointsY
                }
            }
            printf("Time: %.2f percent completed\n", 100.*stepNumber*dt/T);
	        counterOutput++;
        }
        
        tCurrent += dt;
        stepNumber += 1;

        // swapping time-layers
        // real *tmp;

        ////// swap V
        tmp = VOld; VOld = VNew; VNew = tmp;
        ///// swap m
        tmp = mOld; mOld = mNew; mNew = tmp;
        // ///// swap n
        //tmp = nOld; nOld = nNew; nNew = tmp;
        //// swap h
        tmp = hOld; hOld = hNew; hNew = tmp;
        ///// swap p
        tmp = pOld; pOld = pNew; pNew = tmp;
        ///// swap q
        tmp = qOld; qOld = qNew; qNew = tmp;
        ///// swap d
        tmp = dOld; dOld = dNew; dNew = tmp;
        ///// swap f
        tmp = fOld; fOld = fNew; fNew = tmp;
        /*
        ///// swap INa
        tmp = INaOld; INaOld = INaNew; INaNew = tmp;
        ///// swap IK
        tmp = IKOld; IKOld = IKNew; IKNew = tmp;
        ///// swap ILeak
        tmp = ILeakOld; ILeakOld = ILeakNew; ILeakNew = tmp;
        ///// swap IHyperpolar
        tmp = IHyperpolarOld; IHyperpolarOld = IHyperpolarNew; IHyperpolarNew = tmp;
        ///// swap ISlow
        tmp = ISlowOld; ISlowOld = ISlowNew; ISlowNew = tmp;
        */
        /* ((stepNumber % (10* 5000)) == 0) {
            #pragma acc update host(VOld[0:numPointsTotal]) 
	        for(const auto& variable: variables) {
                Write2VTK(variable.first, numPointsX, variable.second, hx, counterOutput); // for now: numPointsX == numPointsY
                //Write2VTK("V", numPointsX, variables["V"], hx, counterOutput); // for now: numPointsX == numPointsY
            }
            printf("Step #%d is performed\n", stepNumber);
	        counterOutput++;
        }*/


    } // while (tCurrent < T)


} // acc data

    printf("\nCalculations finished. Elapsed time = %.2e sec\n", ((real)(clock() - start))/CLOCKS_PER_SEC);

    
    // cleaning up
    delete[] VOld;
    delete[] VNew;
    delete[] mOld;
    delete[] mNew;
    //delete[] nOld;
    //delete[] nNew;
    delete[] hOld;
    delete[] hNew;
    delete[] pOld;
    delete[] pNew;
    delete[] qOld;
    delete[] qNew;
    delete[] dOld;
    delete[] dNew;
    delete[] fOld;
    delete[] fNew;


    return 0;
}
