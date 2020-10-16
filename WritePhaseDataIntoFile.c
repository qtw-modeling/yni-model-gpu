#include "state.h"
#include "common.h" // one for all includes
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
//#define T (2000.) //(1000.) // old val: 500 // endtime (in ms)
#define dt 0.005 // old val = 1e-4 // timestep (in ms)

// model parameters
#define Cm 1.
#define VRest (-60.) // NOTE: there exists no resting potential for SA node

// tissue parameters
//#define Dx (7e-3) //88e-3 //1e-3 //1e-3 // conductivity
//#define Dy (7e-3) //88e-3 //1e-3 //1e-3 // conductivity

// Currents are in the end of enum
enum vars {V_, m_, h_, p_, q_, d_, f_, INa_, IK_, ILeak_, IHyperpolar_, ISlow_};


real* MemAlloc(int n)
{
    return (real*)malloc(n * sizeof(real));
}


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
int main()
{
    // we dont need to perform memalloc for members within these structs
    struct State* Old = (struct State*)malloc(sizeof(struct State));
    struct State* New = (struct State*)malloc(sizeof(struct State));
    
    // initial conditions: only for Old structs: New will be calculated in the loop
    Old->V = VRest;
    Old->m = 0.067;                //m_inf_CPU(VRest); // 0.5
    Old->h = 0.999; //h_inf_CPU(VRest); // 0.5
    Old->p = 0.;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->q = 0.;    // may be false; TODO perform calcs with higher T
    Old->d = 0.; //0.;
    Old->f = 1.; //1.;

    real* stateOfPhase = MemAlloc(7); // 7 --- number of vars., wiht currents excluded
    
    // dont forget to memalloc pointers within Currents struct;
    // we use only single cell
    struct Currents* CurrentsOld = (struct Currents*)malloc(sizeof(struct Currents)); //new Currents;
    
    // ??? maybe we dont need CurrentsNew??? we are not swapping OLD and NEW
    struct Currents* CurrentsNew = (struct Currents*)malloc(sizeof(struct Currents));
    // dont forget to memalloc pointers within Currents struct;
    // we use only single cell
    CurrentsOld->INa = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->IK = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->ILeak = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->IHyperpolar = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->ISlow = (real*)malloc(sizeof(real)); //new real[1];

    // ...same for "New" struct here
    CurrentsNew->INa = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->IK = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->ILeak = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->IHyperpolar = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->ISlow = (real*)malloc(sizeof(real)); //new real[1];
    

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

                // remember values, corresponding to (t0): below, we will begin the timestepping from this state
                stateOfPhase[V_] = New->V;
                stateOfPhase[m_] = New->m;
                stateOfPhase[h_] = New->h;
                stateOfPhase[p_] = New->p;
                stateOfPhase[q_] = New->q;
                stateOfPhase[d_] = New->d;
                stateOfPhase[f_] = New->f;
            }

        }

        
        tCurrent += dt;
        counter += 1;

        // swapping time layers
        struct State* tmp;
        tmp = Old; Old = New; New = tmp;

        //printf("Iteration #%d\n", counter);

    } // while

    //printf("t0 = %.2f, t1 = %.2f\n", time0, time1);
    //std::cin.get();

    // set vars, calculated within the loop
    real period = (time1 - time0); //*0.5; // period of oscillations; remove "0.5" when period calc bug is found!
    //printf("First loop is finished; period of oscillations: %.2f ms\n", period);
    
    
    // repeat the loop (calculations) again and find V(phi)
    //tCurrent = 0; // again
    
    //real tOfPhase = TimeViaPhase(phase, period, time0);
    //printf("Phase: %.2f, tOfPhase: %.2f\n", phase, tOfPhase);
    //std::cin.get();

    //real VOfPhase; // to be determined in the loop below
    //std::map<std::string, real> stateOfPhase;

    // initial conditions, corresponding to (t0): only for "Old" structs, "New" ones will be calculated in the loop
    Old->V = stateOfPhase[V_];
    Old->m = stateOfPhase[m_]; // 0.067; //m_inf_CPU(VRest); // 0.5
    Old->h = stateOfPhase[h_]; // 0.999; //h_inf_CPU(VRest); // 0.5
    Old->p = stateOfPhase[p_]; // 0.;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->q = stateOfPhase[q_]; // 0.;    // may be false; TODO perform calcs with higher T
    Old->d = stateOfPhase[d_]; // 0.;    //0.;
    Old->f = stateOfPhase[f_]; // 1.;    //1.;


    int stepNumberAfterTime0 = 0;
    
    //FILE* ff = fopen("phase_data_yni_model.txt", "w"); // 4writing phase data in a file
    unsigned char buffer[10];
    FILE* fp = fopen("phase_data_yni_model.bin", "wb");  // r for read, b for binary
    //fwrite(buffer, sizeof(buffer), 1, ptr); // read 10 bytes to our buffer

    // (again): main loop: timestepping, we begin from (t0)
    tCurrent = time0;
    while (tCurrent <= time1)
    {
        if ( /* (tCurrent >= time0)  && */ (stepNumberAfterTime0 % 10 == 0)) // print to file each 5 msec: 5/dt = 1000
        {    
            real PHASE = 2*M_PI * (tCurrent - time0) / period; // standart PHASE calc formula

            //VOfPhase = Old->V; // nearest-neighbour iterpolation; change to linear!
            stateOfPhase[V_] = Old->V;
            stateOfPhase[m_] = Old->m;
            stateOfPhase[h_] = Old->h;
            stateOfPhase[p_] = Old->p;
            stateOfPhase[q_] = Old->q;
            stateOfPhase[d_] = Old->d;
            stateOfPhase[f_] = Old->f;

            // values, corresponding to the PHASE; writing values to the file in a SINGLE line
            fwrite(&PHASE, sizeof(real), 1, fp); // fprintf(ff, "%.2f,", PHASE);
            // then, write state values
            fwrite(stateOfPhase, 7*sizeof(real), 1, fp) ;
            //fclose(fp);


            //stepNumberAfterTime0 += 1; // this line is INCORRECT! The counter should be updated below (after formulas for making timestep)

            printf("Iteration #%d\n", stepNumberAfterTime0);
            
            //break;
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
        stepNumberAfterTime0 += 1; // UNCOMMENT THIS

        // swapping time layers
        struct State* tmp;
        tmp = Old; Old = New; New = tmp;

    } // while


    fclose(fp);

/*
    // DEBUG: checking reading the binary file
    if( (fp = fopen("phase_data_yni_model.bin", "rb")) == NULL) 
    {
        printf("Cannot open file");
        return 1;
    }

    real PHASE4Debug;
    real* stateOfPhase4Debug = (real*)malloc(7*sizeof(real));
    
for (int ii = 0; ii <= 20; ii++)
{
    // чтение за раз всего массива balance
    //for (int ii = 0; ii < SMTH; ii++)
    //{
        //fread(&PHASE4Debug, sizeof(real), 1, fp);
        fread(stateOfPhase4Debug, 7*sizeof(real), 1, fp);
    //}

    // вывод содержимого массива
    //printf("PHASE = %.2f\n", PHASE4Debug);
    for(int i = 0; i <= 6; i++) 
        printf("%.2f\n", stateOfPhase4Debug[i]);
}
    
    
    // the phase of next chuck
    //fread(&PHASE4Debug, sizeof(real), 1, fp);
    fread(stateOfPhase4Debug, 7*sizeof(real), 1, fp);
    
    //printf("PHASE = %.2f\n", PHASE4Debug);
    for(int i = 0; i <= 6; i++) 
        printf("%.2f\n", stateOfPhase4Debug[i]);

    
    fclose(fp);
*/

    return 0;
}