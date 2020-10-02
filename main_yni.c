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
#define hx 0.07 //1. // uncomment if cells are connected // (1./numSegmentsX)
#define hy 0.07 //1. // uncomment if cells are connected // (1./numSegmentsY)
//#define T (5000.) //(1000.) // old val: 500 // endtime (in ms)
#define dt 0.005 // old val = 1e-4 // timestep (in ms)

// model parameters
#define Cm 1.
#define VRest (-60.) // NOTE: there exists no resting potential for SA node

// tissue parameters
#define Dx (6*7e-3) //88e-3 //1e-3 //1e-3 // conductivity
#define Dy (6*7e-3) //88e-3 //1e-3 //1e-3 // conductivity


// Currents are in the end of the enum
enum vars {V_, m_, h_, p_, q_, d_, f_, INa_, IK_, ILeak_, IHyperpolar_, ISlow_};


real* MemAlloc(int n)
{
    return (real*)malloc(n * sizeof(real));
}


void Write2VTK(const int n, real* p, const real h, const int step)
{
    char fn[256];
    sprintf(fn, "./output/V.%d.vtk", step); 

    // C code
    FILE* f = fopen(fn, "w");
    fprintf(f, "# vtk DataFile Version 3.0\nSolution\nASCII\nDATASET RECTILINEAR_GRID\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", n + 1, n + 1);
    fprintf(f, "X_COORDINATES %d double\n", n + 1);

    int i;
    for (i = 0; i < n +1; i++)
    {
        fprintf(f, "%.2e ", i*h);
    }
    fprintf(f, "\n");

    fprintf(f, "Y_COORDINATES %d double\n", n + 1);
    for (i = 0; i < n +1; i++)
    {
        fprintf(f, "%.2e ", i*h);
    }
    fprintf(f, "\n");

    fprintf(f, "Z_COORDINATES 1 double\n0\n");
    fprintf(f, "CELL_DATA %d\n", (n * n));
    fprintf(f, "SCALARS V_membrane double\nLOOKUP_TABLE default\n");    
    
    int j;
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++)
            fprintf(f, "%.2e ",  p[j * n + i]);
        fprintf(f, "\n");
    }

    fclose(f);


/* C++ code
    std::fstream f(fn, std::ios::out);
    f << "# vtk DataFile Version 3.0" << std::endl;
    f << "Solution" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET RECTILINEAR_GRID" << std::endl;
    f << "DIMENSIONS " << n + 1 << " " << n + 1 << " 1" << std::endl;
    f << "X_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Y_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Z_COORDINATES 1 double\n0" << std::endl;
    f << "CELL_DATA " << (n * n) << std::endl;
    f << "SCALARS V_membrane double\nLOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++)
            f << p[j * n + i] << " ";
        f << std::endl;
    }
    f.close();
    */
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
        

    }
    
    fclose(fp);
    
    return stateOfPhaseFinal; //VOfPhase;
}


void SetInitialConditions_CPU(real* V, real* m, real* h, real* p, real* q, real* d, 
real* f, real value, int numPointsX, int numPointsY) {
    int idx;
    //std::srand(unsigned(1.)); // initial seed for random number generator
    //real randomNumber;

    // single initial peak
    //int iCenter = (int)((real)numSegmentsX /2.);
    //int jCenter = (int)((real)numSegmentsY /2.); // tmp values
    //int idxCenter = CalculateLinearCoordinate_CPU(iCenter, jCenter);

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) {

            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            //randomNumber =  ((real)(std::rand() % 20))/20.; // 4phase setting

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

                real* stateForPhase = CalculateStateFromFile(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase["V"]);
                //std::cin.get();

                V[idxCenter] = stateForPhase[V_];  //VviaPhase(phase); //M_PI/12. //VRest;
                m[idxCenter] = stateForPhase[m_]; //0.067;//m_inf_CPU(VRest); // 0.5
                h[idxCenter] = stateForPhase[h_]; //0.999; //h_inf_CPU(VRest); // 0.5
                p[idxCenter] = stateForPhase[p_]; //0.;// 0.1; // may be false; TODO perform calcs with higher T 
                q[idxCenter] = stateForPhase[q_]; //0.; // may be false; TODO perform calcs with higher T  

                d[idxCenter] = stateForPhase[d_]; //0.; //0.;
                f[idxCenter] = stateForPhase[f_]; //1.; //1.;

                // for progress checking: in percents
                // printf("Set. initial cond: %.2f percent completed\n", 
                //        100.*idxCenter / CalculateLinearCoordinate_CPU(numPointsX - 1, numPointsY - 1, numPointsX));
            }

    // after filling the whole area: "fill" borders wiht Neumann boundary cond.
    // the borders: Neumann boundary conditions
    for (int j = 0; j < numPointsY; j++)
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


    const char* T_s = argv[2];
    const real T = atof(T_s); // atoi(argv[2]);

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


    // initializing before timesteppin'
    SetInitialConditions_CPU(VOld, mOld, /* nOld,*/ hOld, pOld, qOld, dOld, fOld, 0., numPointsX, numPointsY);
    //SetInitialConditions_CPU(VNew, mNew, /* nNew,*/ hNew, pNew, qNew, dNew, fNew, 0.); // for avoiding "junk" values in all '...New' arrays

    real tCurrent = 0.;
    int stepNumber = 0;
    int counterOutput = 1;

    printf("Timesteppin' begins...\n");
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
                
                // inner cells
                //if (i >= 1 && j >= 1 && i <= (numPointsX - 1) && j <= (numPointsY - 1))
                //{
                    // for short names
                    int idxUp = CalculateLinearCoordinate(i, j + 1, numPointsX);
                    int idxDown = CalculateLinearCoordinate(i, j - 1, numPointsX);
                    int idxLeft = CalculateLinearCoordinate(i - 1, j, numPointsX);
                    int idxRight = CalculateLinearCoordinate(i + 1, j, numPointsX);

                    // reading VOld[idxCenter] ONLY ONCE, 4speedup
                    //real VOldCenter = VOld[idxCenter];

                    ///////////////// gating variables: ode ("reaction") step

                    // TODO: make only ONE read of VOld[idxCenter], etc from memory; to more speedup, esp. for GPU
                    mNew[idxCenter] = m_inf(VOld[idxCenter]) + (mOld[idxCenter] - m_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_m(VOld[idxCenter]) + beta_m(VOld[idxCenter])));


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
                    )
                    + dt / Cm * (TotalIonCurrent(idxCenter, VOld[idxCenter], mOld[idxCenter],
                                                             hOld[idxCenter], pOld[idxCenter], 
                                                            qOld[idxCenter], dOld[idxCenter], fOld[idxCenter],
                                                            INaOld, IKOld, ILeakOld, IHyperpolarOld, ISlowOld)
                                                                         + I_Stim(i, j, 0.*1e0)); // "standart" I_stim = 1e0;
                    
               //} // if
               
               

               /*
               // the borders: Neumann boundary conditions
               else
               {
                    int idxNear;
                    
                    if ((i == 0) && (j >= 1) && (j <= numPointsY - 1)) // left border, except for corner cells
                        idxNear = CalculateLinearCoordinate(i + 1, j, numPointsX);
                    else if ((j == 0) && (i >= 1) && (i <= numPointsX - 1)) // bottom, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j + 1, numPointsX);
                    else if ((j == numPointsY) && (i >= 1) && (i <= numPointsX - 1)) // top, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j - 1, numPointsX);
                    else if ((i == numPointsX) && (j >= 1) && (j <= numPointsY - 1)) // right, except for corner cells
                        idxNear = CalculateLinearCoordinate(i - 1, j, numPointsX);
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

            */

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


        if ( (stepNumber % 2000) == 0)  // output each 10 msec: 10/dt = 2000
        {
            //if ( (stepNumber) % (int)(T/dt/500)  == 0 ) {
            #pragma acc update host(VOld[0:numPointsTotal], mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
            pOld[0:numPointsTotal], qOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], \
            INaOld[0:numPointsTotal], IKOld[0:numPointsTotal], ILeakOld[0:numPointsTotal], IHyperpolarOld[0:numPointsTotal], \
            ISlowOld[0:numPointsTotal])
            
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
            for (int k = 0; k < NUM_VARS; k++)
            {    
                Write2VTK_2D_noGhosts(numPointsX, variablesOld[k], hx, outNumber, k); // for now: numPointsX == numPointsY
                
                //Write2VTKWithGhosts(numPointsX, variablesOld[V_], hx, outNumber); // for now: numPointsX == numPointsY
                //Write2VTK_2D_noGhosts(numPointsX, variablesOld[V_], hx, outNumber, k); // for now: numPointsX == numPointsY
                //Write2VTK_2D_noGhosts(numPointsX, variablesOld[INa_], hx, outNumber, "INa"); // for now: numPointsX == numPointsY
            }
            
            //Write2VTK_2D(numPointsX, variables[V_], hx, counterOutput); // for now: numPointsX == numPointsY

            printf("Progress: %.2f %% completed\n", 100.*stepNumber*dt/T);
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

    real elapsedTime = (real)( ((real)(clock() - start))/CLOCKS_PER_SEC );
    printf("\nCalculations finished. Elapsed time: %.2e s.\n", elapsedTime);


    // OUTPUT: elapsed time
    // storing output file's name in char[]
    /* string */ // char tmp_output_file[256];
    // sprintf(tmp_output_file, argv[4]);
    char* outFile4Timing = concat( concat(concat(concat("timing_dim_", numCellsX_s), "x"), numCellsX_s), 
                            concat(concat("_cells_T_", T_s), "_ms.txt"));

    // printing elapsed time into a file
    FILE* ff = fopen(concat("./timings/", outFile4Timing), "w"); // "timings" --- folder name
    fprintf(ff, "%.6e", elapsedTime); // in sec.
    fclose(ff);
    
    // cleaning up
    for (int k = 0; k < NUM_VARS; k++)
    {
        free(variablesOld[k]);
        free(variablesNew[k]);
    }

    /*
    // cleaning up OLD VERS
    free(VOld);
    free(VNew);
    free(mOld);
    free(mNew);
    free(hOld);
    free(hNew);
    free(pOld);
    free(pNew);
    free(qOld);
    free(qNew);
    free(dOld);
    free(dNew);
    free(fOld);
    free(fNew);
    */

    return 0;
}
