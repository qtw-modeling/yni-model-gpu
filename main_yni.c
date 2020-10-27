//
// Created by QuoTheWhite on 27/03/2019.
//
//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <fstream>
//#include <cstdlib>
//#include "openacc.h"


//#include <string>
//#include <map>
//#include <cmath>
#include <math.h> // 4some specific constants

// 4 C++11 includes
//#include <vector>

//typedef double real;
#include "state.h"
#include "common.h" // one for all includes
#include "extra.h"
#include "CurrentNa.h"
#include "CurrentK.h"
#include "CurrentLeak.h"
#include "CurrentHyperpolar.h"
#include "CurrentSlow.h"


// grid parameters
//#define numSegmentsX 30 //10
//#define numSegmentsY 30 //10
//#define numPointsX (numSegmentsX + 1) // = numCells
//#define numPointsY (numSegmentsY + 1) // = numCells
//#define numPointsTotal (numPointsX * numPointsY)
//#define hx 0.07 //1. // uncomment if cells are connected // (1./numSegmentsX)
//#define hy 0.07 //1. // uncomment if cells are connected // (1./numSegmentsY)
//#define T (5000.) //(1000.) // old val: 500 // endtime (in ms)
#define dt 0.005 // 0.005 --- old "main" val // old val = 1e-4 // timestep (in ms)

// model's parameters
#define Cm 1. // in muF/cm^2; is surface density 
#define CM (20e-6) // is abs. value; in muF; original value == 20 pF
#define VRest (-60.) // NOTE: there exists no resting potential for SAN's cells

#define GX (7e-5) // (6e-5) //  in mS; 6e-5 mS --- is equal to 60 nS in article SAN-fibros(2018)
#define GY (7e-5) // (6e-5) // in mS; is equal to 60 nS in article SAN-fibros, 2018

//#define CELL_SIZE (7e-3) // in cm; is equal to 70 mu*m in article SAN-fibros(2018)
#define CELL_DIAM (5e-3) // (7e-3) // in cm; is equal to 70 mu*m in article SAN-fibros(2018); square cell

// tissue parameters
#define DX (GX*CELL_DIAM*CELL_DIAM/CM)  // (60e-3) // (6*7e-3)  // conductivity
#define DY (GY*CELL_DIAM*CELL_DIAM/CM) // (60e-3) // (6*7e-3)  // conductivity

#define AREA_SIZE_X (0.5) // (1.4) // in cm; from article SAN-fibros(2018)
#define AREA_SIZE_Y (0.5) // (1.4) // in cm; from article SAN-fibros(2018)

// Currents are in the end of the enum
enum vars {V_, m_, h_, p_, q_, d_, f_, INa_, IK_, ILeak_, IHyperpolar_, ISlow_};


real* MemAlloc(int n)
{
    return (real*)malloc(n * sizeof(real));
}


int CalculateLinearCoordinate_CPU(int i, int j, int numPointsX) {
    return i + j*numPointsX;
}

#pragma acc routine
int CalculateLinearCoordinate(int i, int j, int numPointsX) {
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
real* CalculateStateFromFile(real phase) 
{
    
    FILE* fp;
    /* открытие для чтения */
    if( (fp = fopen("phase_data_yni_model.bin", "rb")) == NULL) 
    {
        printf("Cannot open file");
        exit(1); //return 1;
    }

    real phaseLocalOld, phaseLocalNew;
    real* stateOfPhaseOld = (real*)malloc(7*sizeof(real));
    real* stateOfPhaseNew = (real*)malloc(7*sizeof(real));
    real* stateOfPhaseFinal = (real*)malloc(7*sizeof(real));
    real* tmp; // 4switching "timelayers"

    // setting "initial conditions"
    fread(&phaseLocalOld, sizeof(real), 1, fp);
    /* reading the whole array in one step */
    fread(stateOfPhaseOld, 7*sizeof(real), 1, fp);
    
    // going over all data in the bin file
    //int counter = 1;
    while (1)
    {
        
        fread(&phaseLocalNew, sizeof(real), 1, fp);
        /* reading the whole array in one step */
        fread(stateOfPhaseNew, 7*sizeof(real), 1, fp);
        
        if (phaseLocalNew >= phase)
        {
            //break; // stands for right-neighbour interpolation, since we have "passed" over "phase"
        
            for (int ii = 0; ii <= 6; ii++)
                stateOfPhaseFinal[ii] = CalculateLinearInterpolate(phase, phaseLocalOld, phaseLocalNew, stateOfPhaseOld[ii], stateOfPhaseNew[ii]);
        
            break;
        }

        tmp = stateOfPhaseOld; stateOfPhaseOld = stateOfPhaseNew; stateOfPhaseNew = tmp;
        
        //counter++;
        

        //printf("Iteration #%d over the binary file\n", counter);

    }
    
    fclose(fp);
    
    return stateOfPhaseFinal; //VOfPhase;
}


void SetInitialConditions_CPU(real* V, real* m, real* h, real* p, real* q, real* d, 
real* f, real value, int numPointsX, int numPointsY, real hx, real hy) 
{
    int idx;

    for (int j = 0; j < numPointsY; j++)
    {
        for (int i = 0; i < numPointsX; i++)
        {

            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            //randomNumber =  ((real)(std::rand() % 20))/20.; // 4phase setting

            // the borders: Dirichlet boundary conditions
            //if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1)) {
                // TODO: find out about the values
                
                // TODO: wrong formulas below ////////////////////////////////////////////////////
                //int jTilde = numPointsY/2 - j; // j = j^{star}; jTilde = jTilde_{star}
                //int iTilde = i - (numPointsX/2); // i = i^{star}; iTilde = iTilde_{star}

                // for phase calculation: using angle in polar coords

                /*
                // case 1: discrete media and nechetniy number of cells
                real LTotal = numPointsX*hx; // should = numPointsY*hy
                real L = (j*hy + hy/2.) - (numPointsY*hy/2.); // LTotal/2. - (j*hx + hx/2.) --- old formula, when j-order was incorrect
                real lsmall = LTotal/2. - ( (numPointsX - 1 - i)*hx + hx/2.);

                real phase = atan2(L, lsmall); // = angle in polar coords; use atan2() func instead of atan() !
                
                //phase = phase - M_PI/2.; // phase shift according to R.Syunyaev article, 2017
                
                // check sign: atan2() returns vals. from [-pi, pi]
                if (phase < 0)
                {
                    phase += 2*M_PI;
                }
                */

                // case 2: continious media and chetniy (as req in GPU-calcs) number of "cells";
                // ghost "cells" ARE INCLUDED in the process
                real LTotal = (numPointsX - 2)*hx; // same should be in Y-axis
                real L = (j*hy + hy/2.) - (LTotal/2. + hy);
                real lsmall = (i*hx + hx/2.) - (LTotal/2. + hx);
                real phase = atan2(L, lsmall); // = angle in polar coords; use atan2() func instead of atan() !
                
                //phase = phase - M_PI/2.; // phase shift according to R.Syunyaev article, 2017
                
                // check sign: atan2() returns vals. from [-pi, pi]
                if (phase < 0)
                {
                    phase += 2*M_PI;
                }

                //printf("Phase = %.2f deg.\n", phase*180/M_PI);
                //std::cin.get();
                // TODO //////////////////////////////////////////////////////////////

                real* stateForPhase = CalculateStateFromFile(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase[V_]);
                //std::cin.get();

                V[idxCenter] = stateForPhase[V_];  //VviaPhase(phase); //M_PI/12. //VRest;
                m[idxCenter] = stateForPhase[m_]; //0.067;//m_inf_CPU(VRest); // 0.5
                h[idxCenter] = stateForPhase[h_]; //0.999; //h_inf_CPU(VRest); // 0.5
                p[idxCenter] = stateForPhase[p_]; //0.;// 0.1; // may be false; TODO perform calcs with higher T 
                q[idxCenter] = stateForPhase[q_]; //0.; // may be false; TODO perform calcs with higher T  

                d[idxCenter] = stateForPhase[d_]; //0.; //0.;
                f[idxCenter] = stateForPhase[f_]; //1.; //1.;

                // for progress checking: in percents
                //printf("Set. initial cond: %.2f percent completed\n", 
                //        100.*idxCenter / CalculateLinearCoordinate_CPU(numPointsX - 1, numPointsY - 1, numPointsX));
            }
    }

    // after filling the whole area: "fill" borders wiht Neumann boundary cond.
    // the borders: Neumann boundary conditions
    for (int j = 0; j < numPointsY; j++)
    {    
        for (int i = 0; i < numPointsX; i++)
        {
            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            
            // borrder cells, including corner cells
            if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsX - 1))
            {
                int idxNear;

                if ((i == 0)) //&& (j >= 1) && (j <= numSegmentsY - 1)) // left border, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i + 1, j, numPointsX);
                if ((j == 0)) //&& (i >= 1) && (i <= numSegmentsX - 1)) // bottom, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j + 1, numPointsX);
                if ((j == numPointsY - 1)) // && (i >= 1) && (i <= numSegmentsX - 1)) // top, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j - 1, numPointsX);
                if ((i == numPointsX - 1)) // && (j >= 1) && (j <= numSegmentsY - 1)) // right, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i - 1, j, numPointsX);

                // what about corner cells? for now, they are not treated (?)
                V[idxCenter] = V[idxNear];
                m[idxCenter] = m[idxNear];
                h[idxCenter] = h[idxNear];
                p[idxCenter] = p[idxNear];
                q[idxCenter] = q[idxNear];
                d[idxCenter] = d[idxNear];
                f[idxCenter] = f[idxNear];
            } // if
        } // for i
    } // for j

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



int main(int argc, char** argv) 
{

    // setting a GPU for the computations; NOTE: req "openacc.h"!
    // acc_set_device_num(1, acc_device_nvidia);

    // we pass number of cells as command line args; 
    // reading the params from the console:
    const char* numCellsX_s = argv[1];
    int numCellsX = atoi(numCellsX_s); // atoi(argv[1]);
    int numCellsY = numCellsX; // atoi(argv[2]);
    int numCellsTotal = numCellsX*numCellsY;
    // int serieOfLaunchesNum = atoi(argv[3]);
    int numPointsX = numCellsX + 2; // includes 2 ghost cells; numCells + 2 = numPointsX
    int numPointsY = numCellsY + 2; // same
    int numPointsTotal = numPointsX * numPointsY;

    real hx = AREA_SIZE_X / numCellsX; // AREA_SIZE_X / (numPointsX - 1); TODO: check
    real hy = AREA_SIZE_Y / numCellsY; // AREA_SIZE_Y / (numPointsY - 1); TODO: check

    const char* T_s = argv[2];
    const int T = atoi(T_s); // atoi(argv[2]);


    const char* mode = argv[3];
    // allocating memory
    
    // C++11 style alloc
    //std::vector<real> vecVOld(numPointsTotal);
    //real* VOld = &vecVOld.front();

    // gating vars
    real* VOld = MemAlloc(numPointsTotal);
    real* mOld = MemAlloc(numPointsTotal);
    real* hOld = MemAlloc(numPointsTotal);
    real* pOld = MemAlloc(numPointsTotal);
    real* qOld = MemAlloc(numPointsTotal);
    real* dOld = MemAlloc(numPointsTotal);
    real* fOld = MemAlloc(numPointsTotal);

    real* VNew = MemAlloc(numPointsTotal);
    real* mNew = MemAlloc(numPointsTotal);
    real* hNew = MemAlloc(numPointsTotal);
    real* pNew = MemAlloc(numPointsTotal);
    real* qNew = MemAlloc(numPointsTotal);
    real* dNew = MemAlloc(numPointsTotal);
    real* fNew = MemAlloc(numPointsTotal);

    // currents
    real* INaOld = MemAlloc(numPointsTotal);
    real* IKOld = MemAlloc(numPointsTotal);
    real* ILeakOld = MemAlloc(numPointsTotal);
    real* IHyperpolarOld = MemAlloc(numPointsTotal);
    real* ISlowOld = MemAlloc(numPointsTotal);

    real* INaNew = MemAlloc(numPointsTotal);
    real* IKNew = MemAlloc(numPointsTotal);
    real* ILeakNew = MemAlloc(numPointsTotal);
    real* IHyperpolarNew = MemAlloc(numPointsTotal);
    real* ISlowNew = MemAlloc(numPointsTotal);

    real* tmp; // a pointer for swapping time-layers 'n' and 'n+1'

    const int NUM_VARS = 12;
    // for output in a loop
    real* variablesOld[NUM_VARS];
    variablesOld[V_] = VOld;
    
    variablesOld[m_] = mOld;
    variablesOld[h_] = hOld;
    variablesOld[p_] = pOld;
    variablesOld[q_] = qOld;
    variablesOld[d_] = dOld;
    variablesOld[f_] = fOld;
    
    variablesOld[INa_] = INaOld;
    variablesOld[IK_] = IKOld;
    variablesOld[ILeak_] = ILeakOld;
    variablesOld[IHyperpolar_] = IHyperpolarOld;
    variablesOld[ISlow_] = ISlowOld;
    
    real* variablesNew[NUM_VARS];
    variablesNew[V_] = VNew;
    
    variablesNew[m_] = mNew;
    variablesNew[h_] = hNew;
    variablesNew[p_] = pNew;
    variablesNew[q_] = qNew;
    variablesNew[d_] = dNew;
    variablesNew[f_] = fNew;
    
    variablesNew[INa_] = INaNew;
    variablesNew[IK_] = IKNew;
    variablesNew[ILeak_] = ILeakNew;
    variablesNew[IHyperpolar_] = IHyperpolarNew;
    variablesNew[ISlow_] = ISlowNew;
    
    //char varNames[NUM_VARS] = {"V", "m", "h", "p", "q", "d", "f"};
    // = {"V", "m", "h", "p", "q", "d", "f"};


    // initializing before timestepping
    SetInitialConditions_CPU(VOld, mOld, /* nOld,*/ hOld, pOld, qOld, dOld, fOld, 0., numPointsX, numPointsY, hx, hy);
    //SetInitialConditions_CPU(VNew, mNew, /* nNew,*/ hNew, pNew, qNew, dNew, fNew, 0.); // for avoiding "junk" values in all '...New' arrays

    // various output params
    real tCurrent = 0.;
    int stepNumber = 0;
    int counterOutput = 1;
    const real timeIntervalOutput = 5.; // in ms
    const int stepsOutput = (int)(timeIntervalOutput/dt); // 1000; // output each 10 ms: 10/dt = 2000; each 5 ms: 5/dt = 1000
    int startOfTimestep; // , timeBetweenOutputs;

    printf("Timestepping begins\n");
    clock_t start = clock();

// pragmas without "-acc" flag --- are ignored?
#pragma acc data copy(VOld[0:numPointsTotal], mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
pOld[0:numPointsTotal], qOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], \
VNew[0:numPointsTotal], mNew[0:numPointsTotal], hNew[0:numPointsTotal], \
pNew[0:numPointsTotal], qNew[0:numPointsTotal], dNew[0:numPointsTotal], fNew[0:numPointsTotal], \
INaOld[0:numPointsTotal], IKOld[0:numPointsTotal], ILeakOld[0:numPointsTotal], IHyperpolarOld[0:numPointsTotal], \
ISlowOld[0:numPointsTotal]), \
deviceptr(tmp)
{
    // main loop: timestepping
    while (tCurrent < T) 
    {
        // for measuring remaining time, which is being printed to console
        if ( ((stepNumber - 1) % stepsOutput == 0) && stepNumber > 0)
            startOfTimestep = clock(); // beginning of measuring time between outputs to .vtk

    // TODO: change order of indexing (i, j)
	#pragma acc parallel \
	present(VOld[0:numPointsTotal], mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
    pOld[0:numPointsTotal], qOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], \
    VNew[0:numPointsTotal], mNew[0:numPointsTotal], hNew[0:numPointsTotal], \
    pNew[0:numPointsTotal], qNew[0:numPointsTotal], dNew[0:numPointsTotal], fNew[0:numPointsTotal]), \
    num_workers(1), vector_length(32)
	{
	
	//#pragma acc loop collapse(2) independent
	#pragma acc loop gang
    for (int j = 1; j < numPointsY - 1; j++)
    {
            #pragma acc loop vector
            for (int i = 1; i < numPointsX - 1; i++) 
            {

                int idxCenter = CalculateLinearCoordinate(i, j, numPointsX);
                
                
                // for short names
                int idxUp = CalculateLinearCoordinate(i, j + 1, numPointsX);
                int idxDown = CalculateLinearCoordinate(i, j - 1, numPointsX);
                int idxLeft = CalculateLinearCoordinate(i - 1, j, numPointsX);
                int idxRight = CalculateLinearCoordinate(i + 1, j, numPointsX);

                // reading VOld[idxCenter] ONLY ONCE, 4speedup
                real VOldCenter = VOld[idxCenter];

                ///////////////// gating variables: ode ("reaction") step

                // TODO: make only ONE read of VOld[idxCenter], etc from memory; to more speedup, esp. for GPU
                mNew[idxCenter] = m_inf(VOldCenter) + (mOld[idxCenter] - m_inf(VOldCenter))
                                                                * exp(-dt * (alpha_m(VOldCenter) + beta_m(VOldCenter)));


                hNew[idxCenter] = h_inf(VOldCenter) + (hOld[idxCenter] - h_inf(VOldCenter))
                                                                * exp(-dt * (alpha_h(VOldCenter) + beta_h(VOldCenter)));
                    
                pNew[idxCenter] = p_inf(VOldCenter) + (pOld[idxCenter] - p_inf(VOldCenter))
                                                                * exp(-dt * (alpha_p(VOldCenter) + beta_p(VOldCenter)));
                   
                qNew[idxCenter] = q_inf(VOldCenter) + (qOld[idxCenter] - q_inf(VOldCenter))
                                                                * exp(-dt * (alpha_q(VOldCenter) + beta_q(VOldCenter)));
                    
                dNew[idxCenter] = d_inf(VOldCenter) + (dOld[idxCenter] - d_inf(VOldCenter))
                                                                * exp(-dt * (alpha_d(VOldCenter) + beta_d(VOldCenter)));
                    
                fNew[idxCenter] = f_inf(VOldCenter) + (fOld[idxCenter] - f_inf(VOldCenter))
                                                                * exp(-dt * (alpha_f(VOldCenter) + beta_f(VOldCenter)));
                   
                    
                //////////////////
                // Euler method's step; no operator splitting is used
                VNew[idxCenter] = VOldCenter
                // uncomment if cells are connected; otherwize --- Nans
                     + dt * (
                            DX / (hx * hx)  * (VOld[idxRight] - 2 * VOldCenter + VOld[idxLeft])
                            + DY /(hy * hy)  * (VOld[idxUp] - 2 * VOldCenter + VOld[idxDown])
                    )
                    + dt / Cm * ( TotalIonCurrent(idxCenter, VOldCenter, mOld[idxCenter],
                                                             hOld[idxCenter], pOld[idxCenter], 
                                                            qOld[idxCenter], dOld[idxCenter], fOld[idxCenter],
                                                            INaOld, IKOld, ILeakOld, IHyperpolarOld, ISlowOld));
                                                                        //  /+ I_Stim(i, j, 0.*1e0) ); // "standart" I_stim = 1e0;
                    
                    // dt / Cm --- is correct! Capacity density (Cm) has to be here
            } // for i
    } // for j
	} // acc parallel 4 inner cells

    #pragma acc parallel async(0) \
    present(VNew[0:numPointsTotal])
    {
        #pragma acc loop vector // seq
        for (int j = 1; j < numPointsY - 1; j++)
        {
            int idxCenterLeftBord = CalculateLinearCoordinate(0, j, numPointsX);
            int idxCenterRightBord = CalculateLinearCoordinate(numPointsX - 1, j, numPointsX);
                            
            int idxNearLeft = CalculateLinearCoordinate(1, j, numPointsX);
            int idxNearRight = CalculateLinearCoordinate(numPointsX - 1 - 1, j, numPointsX);

            VNew[idxCenterLeftBord] = VNew[idxNearLeft];
            VNew[idxCenterRightBord] = VNew[idxNearRight];
        }
    }

    #pragma acc parallel async(1) \
    present(VNew[0:numPointsTotal])
    {
        #pragma acc loop vector // seq
        for (int i = 1; i < numPointsX - 1; i++)
        {
            int idxCenterBottomBord = CalculateLinearCoordinate(i, 0, numPointsX);
            int idxCenterTopBord = CalculateLinearCoordinate(i, numPointsY - 1, numPointsX);
                        
            int idxNearBottom = CalculateLinearCoordinate(i, 1, numPointsX);
            int idxNearTop = CalculateLinearCoordinate(i, numPointsY - 1 - 1, numPointsX);

            VNew[idxCenterBottomBord] = VNew[idxNearBottom];
            VNew[idxCenterTopBord] = VNew[idxNearTop];
        }
    }
                
    #pragma acc wait(0, 1)      


    if ( (stepNumber % stepsOutput /* 1000 */ /* 2000 */) == 0)  // output each 10 ms: 10/dt = 2000; each 5 ms: 5/dt = 1000
    {
        //if ( (stepNumber) % (int)(T/dt/500)  == 0 ) {
        #pragma acc update host(VOld[0:numPointsTotal]) /*, mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
        pOld[0:numPointsTotal], qOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], \
        INaOld[0:numPointsTotal], IKOld[0:numPointsTotal], ILeakOld[0:numPointsTotal], IHyperpolarOld[0:numPointsTotal], \
        ISlowOld[0:numPointsTotal]) */
            
        //variablesOld[V_] = VOld;
        //variablesOld[m_] = mOld;
        //variablesOld[h_] = hOld;
        //variablesOld[p_] = pOld;
        //variablesOld[q_] = qOld;
        //variablesOld[d_] = dOld;
        //variablesOld[f_] = fOld;
        //variablesOld[INa_] = INaOld;
        //variablesOld[IK_] = IKOld;
        //variablesOld[ILeak_] = ILeakOld;
        //variablesOld[IHyperpolar_] = IHyperpolarOld;
        //variablesOld[ISlow_] = ISlowOld;

        int outNumber = stepNumber;
                
        // output 2VTKs in a loop over enums
        //for (int k = 0; k < NUM_VARS; k++)
        //{    
            Write2VTK_2D_noGhosts(numPointsX, variablesOld[V_], hx, outNumber, V_); // for now: numPointsX == numPointsY
            //Write2VTK_2D(numPointsX, variablesOld[V_], hx, outNumber, V_); // for now: numPointsX == numPointsY
        //}

        real timeBetweenOutputs = (real)(clock() - startOfTimestep) / CLOCKS_PER_SEC / 3600. * 60.; // in minutes
        real time4dtAvg = timeBetweenOutputs / (real)stepsOutput; // average only for timesteps between 2 conseq. outputs;
        // != average for the whole "while(t < T)" loop
        
        int numdtRemaining = (int)T/dt - stepNumber;
        printf("Progress: %.2f %% completed || time remaining: %.0f min\n", 100.*stepNumber*dt/T, time4dtAvg * numdtRemaining);
	    counterOutput += 1;
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


        // writing elapsed time
        //real timeOfTimestep = (real)(clock() - startOfTimestep)/CLOCKS_PER_SEC/3600.;
        

    } // while (tCurrent < T)


} // acc data

    real elapsedTime = (real)( ((real)(clock() - start))/CLOCKS_PER_SEC );
    printf("\nCalculations finished. Elapsed time: %.2e s.\n", elapsedTime);


    // OUTPUT: elapsed time
    // storing output file's name in char[]
    /* string */ // char tmp_output_file[256];
    // sprintf(tmp_output_file, argv[4]);
    char* outFile4Timing = concat( concat(concat(concat(concat(mode, "_timing_dim_"), numCellsX_s), "x"), numCellsX_s), 
                            concat(concat("_cells_T_", T_s), "_ms.txt"));

    // printing elapsed time into a file
    FILE* ff;
    if ((ff = fopen(concat("./timings/", outFile4Timing), "w")) == NULL) // "timings" --- folder name
    {
        printf("Cannot open ./timings/-dir to writing. Aborting.");
        exit(1); //return 1;
    }

    fprintf(ff, "%.6e", elapsedTime); // in sec.
    fclose(ff);
    
    /*
    // cleaning up
    for (int k = 0; k < NUM_VARS; k++)
    {
        free(variablesOld[k]);
        free(variablesNew[k]);
    }
    */

    return 0;
}
