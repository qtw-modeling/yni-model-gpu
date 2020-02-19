//
// Created by QuoTheWhite on 27/03/2019.
//
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <cstdlib>
//#include "openacc.h"

// 4C++11 includes
#include <vector>

typedef double real;

// grid parameters
#define numSegmentsX 0
#define numSegmentsY 0
#define numPointsX (numSegmentsX + 1)
#define numPointsY (numSegmentsY + 1)
#define numPointsTotal (numPointsX * numPointsY)
#define hx 1. // uncomment if cells are connected // (1./numSegmentsX)
#define hy 1. // uncomment if cells are connected // (1./numSegmentsY)
#define T 500. // old val: 500 // endtime
#define dt 1e-4 // timestep

// model parameters
#define Cm 1.
#define VRest (-65.)

// tissue parameters
#define Dx 1e-3
#define Dy 1e-3



void Write2VTK(const int n, real* p, const real h, const int step)
{
    char fn[256];
    sprintf(fn, "./output/yni.%d.vtk", step);

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



void SetInitialConditions_CPU(real* V, real* m, real* n, real* h, real value) {
    int idx;
    std::srand(unsigned(1.)); // initial seed for random number generator
    real randomNumber;

    // single initial peak
    int iCenter = (int)((real)numSegmentsX /2.);
    int jCenter = (int)((real)numSegmentsY /2.); // tmp values
    int idxCenter = CalculateLinearCoordinate_CPU(iCenter, jCenter);

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) {

            idx = CalculateLinearCoordinate_CPU(i, j);
            randomNumber =  ((real)(std::rand() % 20))/20.;

            // the borders: Dirichlet boundary conditions
            if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1)) {
                // TODO: find out about the values
                V[idx] = VRest;
                m[idx] = m_inf_CPU(VRest);
                n[idx] = n_inf_CPU(VRest);
                h[idx] = h_inf_CPU(VRest);
            }
            /*else if (idx == idxCenter) { // initial peak
                // TODO: find out about the values
                V[idx] = 0.95; //randomNumber;
                m[idx] = 0.5; //randomNumber;
                n[idx] = 0.;
                h[idx] = 0.;

            }*/
            else {
                //randomNumber =  ((real)(std::rand() % 20))/20.;
                // TODO: find out about the values
                V[idx] = VRest;
                m[idx] = m_inf_CPU(VRest);
                n[idx] = n_inf_CPU(VRest);
                h[idx] = h_inf_CPU(VRest);
            }

        }
}


#pragma acc routine
real I_Stim(int i, int j, real value) {
    //int x0 = (int)((real)numSegmentsX /2.); int y0 = (int)((real)numSegmentsY /2.); // tmp values

    /* uncomment if cells are connected
    if ((i > 1 && i < numPointsX - 2) && (j == 2))
        return value;
    else
        return 0.;
    */

   return value;
}


// ion currents
///////////////////////////////////////////////////
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

#pragma acc routine
real TotalIonCurrent(real V, real m, real n, real h) {
    real iNa = CurrentNa(V, m, n, h);
    real iK = CurrentK(V, m, n, h);
    real iLeak = CurrentLeak(V, m, n, h);

    // TODO: check the sign of the expression below
    return  -(iNa + iK + iLeak);
}




int main() {

    // setting a GPU for the computations; NOTE: req "openacc.h"!
    //acc_set_device_num(1, acc_device_nvidia);

    // allocating memory
    
    // C++11 style alloc
    //std::vector<real> vecVOld(numPointsTotal);
    //real* VOld = &vecVOld.front();
    
    
    real* VOld = new real[numPointsTotal];
    real* mOld = new real[numPointsTotal];
    real* nOld = new real[numPointsTotal];
    real* hOld = new real[numPointsTotal];

    real* VNew = new real[numPointsTotal];
    real* mNew = new real[numPointsTotal];
    real* nNew = new real[numPointsTotal];
    real* hNew = new real[numPointsTotal];

    real* tmp; // a pointer for swapping time-layers 'n' and 'n+1'


    // initializing before timesteppin'
    SetInitialConditions_CPU(VOld, mOld, nOld, hOld, 0.);
    SetInitialConditions_CPU(VNew, mNew, nNew, hNew, 0.); // for avoiding "junk" values in all '...New' arrays

    real tCurrent = 0.;
    int stepNumber = 0;
    int counterOutput = 1;

    printf("Timesteppin' begins...\n");
    clock_t start = clock();


#pragma acc data copy(VOld[0:numPointsTotal], mOld[0:numPointsTotal], nOld[0:numPointsTotal], hOld[0:numPointsTotal], \
		      VNew[0:numPointsTotal], mNew[0:numPointsTotal], nNew[0:numPointsTotal], hNew[0:numPointsTotal]) \
		      deviceptr(tmp)
{
    // main loop: timestepping
    while (tCurrent < T) {

        // TODO: change order of indexing (i, j)
        
	#pragma acc kernels \
	present(VOld[0:numPointsTotal], mOld[0:numPointsTotal], nOld[0:numPointsTotal], hOld[0:numPointsTotal], \
                VNew[0:numPointsTotal], mNew[0:numPointsTotal], nNew[0:numPointsTotal], hNew[0:numPointsTotal])
	{
	
	#pragma acc loop collapse(2) independent
	for (int j = 0; j < numPointsY; j++)
            for (int i = 0; i < numPointsX; i++) {

                int idxCenter = CalculateLinearCoordinate(i, j);
                
                /*
                // uncommnent
                // the borders: Dirichlet boundary conditions
                if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1)) {
                    VNew[idxCenter] = VOld[idxCenter];
                    mNew[idxCenter] = mOld[idxCenter];
                    nNew[idxCenter] = nOld[idxCenter];
                    hNew[idxCenter] = hOld[idxCenter];
                }
                */

                /* else {
                    // for short names
                    int idxUp = CalculateLinearCoordinate(i, j + 1);
                    int idxDown = CalculateLinearCoordinate(i, j - 1);
                    int idxLeft = CalculateLinearCoordinate(i - 1, j);
                    int idxRight = CalculateLinearCoordinate(i + 1, j);

                    */
                    ///////////////// gating variables: ode ("reaction") step
                    mNew[idxCenter] = m_inf(VOld[idxCenter]) + (mOld[idxCenter] - m_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_m(VOld[idxCenter]) + beta_m(VOld[idxCenter])));

                    nNew[idxCenter] = n_inf(VOld[idxCenter]) + (nOld[idxCenter] - n_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_n(VOld[idxCenter]) + beta_n(VOld[idxCenter])));

                    hNew[idxCenter] = h_inf(VOld[idxCenter]) + (hOld[idxCenter] - h_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_h(VOld[idxCenter]) + beta_h(VOld[idxCenter])));


                    //////////////////
                    // "discrete diffusion" step
                    VNew[idxCenter] = VOld[idxCenter]; // 
                    /* uncomment if cells are connected; otherwize --- Nans */ 
                    /* + dt / Cm * (
                            Dx / (hx*hx) * (VOld[idxRight] - 2 * VOld[idxCenter] + VOld[idxLeft])
                            + Dy / (hy*hy) * (VOld[idxUp] - 2 * VOld[idxCenter] + VOld[idxDown])
                    ); */
                    // reaction step
                    VNew[idxCenter] += dt / Cm * (TotalIonCurrent(VOld[idxCenter], mOld[idxCenter],
                                                            nOld[idxCenter], hOld[idxCenter])
                                                                        + I_Stim(i, j, 1e1));

               // } // else
            } // for
	
	} // acc kernels

        tCurrent += dt;
        stepNumber += 1;

        // swapping time-layers
        // real *tmp;

        ////// swap V
        tmp = VOld; VOld = VNew; VNew = tmp;
        ///// swap m
        tmp = mOld; mOld = mNew; mNew = tmp;
        ///// swap n
        tmp = nOld; nOld = nNew; nNew = tmp;
        //// swap h
        tmp = hOld; hOld = hNew; hNew = tmp;


        if ((stepNumber % 5000) == 0) {
            #pragma acc update host(VOld[0:numPointsTotal]) 
	    Write2VTK(numPointsX, VOld, hx, counterOutput); // for now: numPointsX == numPointsY
            printf("Step #%d is performed\n", stepNumber);
	    counterOutput++;
        }


    } // while (tCurrent < T)


} // acc data

    printf("\nCalculations finished. Elapsed time = %.2e sec\n", ((real)(clock() - start))/CLOCKS_PER_SEC);

    // cleaning up
    //delete[] VOld;
    
    delete[] VNew;
    delete[] mOld;
    delete[] mNew;
    delete[] nOld;
    delete[] nNew;
    delete[] hOld;
    delete[] hNew;

    return 0;
}
