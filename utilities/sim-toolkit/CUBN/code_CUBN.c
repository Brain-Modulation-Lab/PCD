/* 
 * @brief mex code to simulate a recurrent random network with excitatory 
 *        and inhibitory LIF CURRENT-BASED neurons. 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
// Needed for building mex file
#include <mex.h>
#include <matrix.h>

#include "ran1.h"  
/* Random number generator routine. To generate real random numbers 0.0-1.0
 * Should be seeded with a negative integer*/
#include "gasdev.h" 
/* gaussian variate generator routine. To generate random number along a normal distribution*/


// Conventions for variable names:
// -------------------------------
// - prefix "e" stands for "excitatory"
// - prefix "i" stands for "inhibitory"
// - prefix "x" stands for "external"
//
// - prefix "a" stands for "AMPA"
// - prefix "g" stands for "GABA"
//
// - prefix "N"   stands for "number of..."
// - prefix "tot" stands for "total"
//
// - prefix "e2e" stands for "excitatory to excitatory"
// - prefix "e2i" stands for "excitatory to inhibitory"
// - prefix "x2e" stands for "external to excitatory"
// - prefix "e2i" stands for "external to inhibitory"
//
// - "T"  stands for greek letter "tau"
//
// - "nrn"  stands for "neurons"
// - "FR"   stands for "firing-rate"
// - "SP"   stands for "spike"
// - "SPC"  stands for "spike-count"
//
// - "rec" stands fo recurrent
// - "ext" stands fo external

/**
 * Main entry called from Matlab
 * @param nlhs number of left-hand-side arguments
 * @param plhs pointers to mxArrays in the left-hand-side
 * @param nrhs number of right-hand-side arguments
 * @param prhs pointers to mxArrays in the right-hand-side
 *
 * This is the entry point of the function
 */

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // Variables declaration ----------------------------------------------
    
    //Time:
    double  Dt;             /* time resolution [in ms]*/
    mwIndex t;              /* time index*/
    double  simulLen;       /* number of steps in the simulation*/
    double  simulLen2;      /* variable to check the length of the second input array, x2iFR*/
    
    mxArray *tmp;
    
    /* Cycles:
     *  -------
     *  Cycles are used to reduce the amount of memory required for storing
     *  spike arrival times. A spike reaches the post-synaptic neuron with
     *  a delay Tl from its time of emission, thus storing arrivals times
     *  would require an array of length simLen for each neuron in the net.
     *  However, in order to save memory, it is possible to use an array of
     *  size Tl/Dt+1 and cycle on the Tl/Dt values.
     * 
     *  EXAMPLE: Let's consider the case Dt=1 and Tl=3. A spike fired at
     *  time t=2 would reach the post-synaptic neuron at time t=5. However,
     *  once the current time is known, we just need to know the arrival
     *  cycle value, cycle=1, to be able to deliver the spike at the right
     *  time.
     * 
     *  t:       0   1   2   3   4   5   6   7   8   9   ...
     *                   |           |
     *                 spike       spike
     *                emission    arrival
     *                   |           |
     *  cycle:   0   1   2   3   0   1   2   3   0   1   ...
     *                   |___________|
     *                         Tl
     * 
     *  Thus a vector of size 4 would be enough for storing all arrival
     *  times of the spikes emitted by one neuron.*/
    
    
    mwSize  eCycSize; /* size of cycle for excitatory nerons*/
    mwSize  iCycSize; /* size of cycle for inhibitory nerons*/
    mwIndex eCycIndx; /* cycle index for excitatory neurons*/
    mwIndex iCycIndx; /* cycle index for inhibitory neurons*/
    mwIndex eCycSP;   /* arrival cycle of spike emitted by pre-syn. ex. nrn*/
    mwIndex iCycSP;   /* arrival cycle of spike emitted by pre-syn. in. nrn*/
    
    /* Number of neurons:*/
    mwSize  eNnrn;   /* number of excitatory neurons*/
    mwSize  iNnrn;   /* number of inhibitory neurons*/
    mwSize  totNnrn; /* total number of neurons*/
    mwIndex nrn;     /* neuron index*/
    mwIndex nrn2;    /*    "     "  */
    
    // Random numbers generation
    double RAND_MAX_double; /* max value returned by function random*/
    double rndNum;          /* stores i.i.d. 0 to 1 distributed random numbers*/
    int seed1;              /* first random seed (it defines connections)*/
    int seed2;              /* second random seed (it defines poisson variate)*/
    
    // Leak membrane potential, V_leaky:
    double V_leaky;
    
    // Spike threshold potential, V_threshold:
    double Vthr;
    
    // Reset potential, V_reset:
    double eVres;
    double iVres;
    
    // Membrane time constants, tau_m:
    double eTm;
    double iTm;
        
    // Refractory period:
    double eTrp;
    double iTrp;
    
    // Synaptic latency, tau_l:
    double eTl;
    double iTl;
    
    // Syaptic currents rise and decay times, tau_r and tau_d:
    double e2eTr;
    double e2eTd;
    double e2iTr;
    double e2iTd;
    double iTr;
    double iTd;
    
    // Synaptic efficacies, J_syn:
    double i2iJ;
    double e2iJ;
    double x2iJ;
    double i2eJ;
    double e2eJ;
    double x2eJ;
    
    /* Membrane resistances, 1/g_leak:*/
    double eRm;
    double iRm;
    
    // Connections:
    double   p;
    mwSize  *A;             // matrix storing neurons connections
    mwSize   a;             // element of A
    mwIndex  con;           // connection index
    mwSize  *Ncon;          // number of outgoing connections of each neuron
    
    // Spikes
    double *tLastSP;        // time last spike was fired 
    
    /*/ Firing rates*/
    double *x2eFR;     /*/ input firing rate to excitatory neuron*/
    double *x2iFR;     /*/ input firing rate to inhibitory neuron*/
    double  exp_x2eFR; /*/ exp(xFR)*/
    double  exp_x2iFR;
    int eCounter;
    int iCounter;
    
    // Membrane potential V
    double *V;
    // Auxiliary variables to compute membrane potential
    double  V_k1;
    double  V_rk;
    double  V_k2;
    double  synI;
    double  synI_rk;
    
    // Synaptic currents I_AMPA, I_GABA
    double  aI; // I_AMPA = I_AMPArec+I_AMPAext
    double  gI;
    
    // Time course of the synaptic currents, s_syn(t)
    double *aS_rec; // recurrent AMPA
    double *aS_ext; // external AMPA
    double *gS;     // GABA
    
    // Auxiliary variables to compute the time course of the synaptic currents
    double C;
    
    double *aX_rec;
    double *aX_ext;
    double  aX_k1;
    double  aX_rk;
    double  aX_k2;
    
    double *gX;
    double  gX_k1;
    double  gX_rk;
    double  gX_k2;
    
    double  aS_k1;
    double  aS_k2;
    double  aS_rk_rec;
    double  aS_rk_ext;
    
    double  gS_k1;
    double  gS_k2;
    double  gS_rk;
    //--------------------------------------------------------------------------
    
    /* Spike counts*/
    double *eSPC;
    double *iSPC;
    double  xSPC;
    double  e2eSPC;
    double  e2iSPC;
    double  i2eSPC;
    double  i2iSPC;

    /*Connections*/
    mwSize *syntemp;
    mwIndex kn, ip;
    int avNcon_outgo;
    int sdNcon_outgo;

    /*OUTPUT variables*/
    double *e2eI;   /*sum of the AMPA currents on excitatory neurons as a function of the time*/
    double *i2eI;   /*sum of the GABA currents on excitatory neurons as a function of the time*/
    unsigned short *eFR;     /* number of excitatory action potentials in each time step*/
    unsigned short *iFR;     /* number of inhibitory action potentials in each time step*/
    
    // end of the variables declaration--------------------------------------------
    
   
    /* Defining the network parameters ------------------------------------
     *
     * Number of neurons:*/
    tmp   = mxGetField(prhs[0], 0, "eNnrn");
    eNnrn = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iNnrn");
    iNnrn = *mxGetPr(tmp);
    
    totNnrn = eNnrn + iNnrn;
    
    // Connection probability:
    tmp   = mxGetField(prhs[0], 0, "p"    );
    p     = *mxGetPr(tmp);
    
    // Membrane time constants [in ms]:
    tmp   = mxGetField(prhs[0], 0, "eTm"  );
    eTm   = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iTm"  );
    iTm   = *mxGetPr(tmp);
    
    // Threshold potential, in [mV]:
    tmp   = mxGetField(prhs[0], 0, "Vthr" );
    Vthr  = *mxGetPr(tmp);
    
    // Reset potential, in [mV]:
    tmp   = mxGetField(prhs[0], 0, "eVres");
    eVres = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iVres");
    iVres = *mxGetPr(tmp);
    
    // Refractory period, in [ms]:
    tmp   = mxGetField(prhs[0], 0, "eTrp" );
    eTrp  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iTrp" );
    iTrp  = *mxGetPr(tmp);
    
    /* Synaptic latency, in [ms]:*/
    tmp   = mxGetField(prhs[0], 0, "eTl"  );
    eTl   = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iTl"  );
    iTl   = *mxGetPr(tmp);
    
    /* Synaptic currents rise and decay times, in [ms]:*/
    tmp   = mxGetField(prhs[0], 0, "e2eTr");
    e2eTr = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "e2eTd");
    e2eTd = *mxGetPr(tmp);
    
    
    tmp   = mxGetField(prhs[0], 0, "e2iTr");
    e2iTr = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "e2iTd");
    e2iTd = *mxGetPr(tmp);
    
    
    tmp   = mxGetField(prhs[0], 0, "iTr"  );
    iTr   = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iTd"  );
    iTd   = *mxGetPr(tmp);
    
    // Synaptic efficacies, in [pA]:
    tmp   = mxGetField(prhs[0], 0, "i2iJ" );
    i2iJ  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "e2iJ" );
    e2iJ  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "x2iJ" );
    x2iJ  = *mxGetPr(tmp);
   
    tmp   = mxGetField(prhs[0], 0, "i2eJ" );
    i2eJ  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "e2eJ" );
    e2eJ  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "x2eJ" );
    x2eJ  = *mxGetPr(tmp);
    
    /* Membrane resistances, in [GOhm]*/
    tmp   = mxGetField(prhs[0], 0, "eRm" );
    eRm  = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iRm" );
    iRm  = *mxGetPr(tmp);
    
    /* Leak membrane potential, in [mV]*/
    tmp   = mxGetField(prhs[0], 0, "V_leaky" );
    V_leaky  = *mxGetPr(tmp);
    
    // end of the definition of the network parameters----------------------
    
    // Printing the network's parameters
    mexPrintf("V_leaky = %.0f mV\n", V_leaky);
    mexPrintf("V_threshold = %.0f mV\n", Vthr);
    mexPrintf("eV_reset = %.0f mV\n", eVres);
    mexPrintf("iV_reset = %.0f mV\n\n", iVres);
    mexPrintf("iT_rise = %.2f ms\n", iTr);    
    mexPrintf("iT_decay = %.2f ms\n", iTd);    
    mexPrintf("e2eT_rise = %.2f ms\n", e2eTr);    
    mexPrintf("e2eT_decay = %.2f ms\n", e2eTd);    
    mexPrintf("e2iT_rise = %.2f ms\n", e2iTr);    
    mexPrintf("e2iT_decay = %.2f ms\n", e2iTd);    
    mexPrintf("eT_latency = %.1f ms\n", eTl);    
    mexPrintf("iT_latency = %.1f ms\n", iTl);    
    mexPrintf("eTm = %.1f ms\n", eTm);    
    mexPrintf("iTm = %.1f ms\n", iTm);   
    mexPrintf("eTrp = %.1f ms\n", eTrp);   
    mexPrintf("iTrp = %.1f ms\n\n", iTrp);    
    mexPrintf("gi2i = %.2f pA\n", i2iJ);    
    mexPrintf("e2eJ = %.2f pA\n", e2eJ);
    mexPrintf("i2eJ = %.2f pA\n", i2eJ);    
    mexPrintf("e2iJ = %.2f pA\n", e2iJ);    
    mexPrintf("x2iJ = %.2f pA\n", x2iJ);   
    mexPrintf("x2eJ = %.2f pA\n\n", x2eJ);
    

    // Integration step, in [ms]:
    tmp = mxGetField(prhs[0], 0, "Dt");
    Dt = *mxGetPr(tmp);
    
    
    // Check of the number of inputs
    if(nrhs<4)
        mexErrMsgTxt("Not enough input arguments.");
    if (nrhs != 5) {
        mexPrintf("Input arguments missing! \n");
    }
   
    // Check of the length of the two external input arrays
    simulLen = mxGetNumberOfElements(prhs[1]); 
    simulLen2 =  mxGetNumberOfElements(prhs[2]);
    if (simulLen != simulLen2) {
        mexErrMsgTxt("excitatory and inhibitory inputs must have the same length \n");    
    }
    // External input arrays
    x2eFR = mxGetPr(prhs[1]);
    x2iFR = mxGetPr(prhs[2]);
        
    /* On the use of the random seeds:
     * -------------------------------
     * The two random seeds are used in the following way:
     * - the first random seed, seed1, affects the network configuration.
     * - the second random seed, seed2, affects the generation of Poisson
     *   random variates, which determine the input to each neuron
     * When seed1=0 or no seed2 is passed as input, time(NULL) is used instead to generate seed1 and seed2.*/
     
     /* Checking if the first random seed has been passed as an input. Using
     * time(NULL) otherwise.*/

    seed1 = *mxGetPr(prhs[3]);
    if(seed1==0)
        seed1 = time(NULL);
    
    // Checking if the second random seed has been passed as an input.
    // Using time(NULL) otherwise.
    if(nrhs==5)
        seed2 = (int) *mxGetPr(prhs[4]);
    else
        seed2 = time(NULL);
    
    if (seed1<0 || seed2 <0){
        mexErrMsgTxt("seed1 and seed2 must be positive integer\n");
    }
   
    /* Normalizing in Dt units --------------------------------------------*/
    eTm   = eTm   / Dt;
    iTm   = iTm   / Dt;
    e2eTr = e2eTr / Dt;
    e2iTr = e2iTr / Dt;
    e2eTd = e2eTd / Dt;
    e2iTd = e2iTd / Dt;
    iTr   = iTr   / Dt;
    iTd   = iTd   / Dt;
    eTrp  = eTrp  / Dt;
    iTrp  = iTrp  / Dt;    
    eTl   = eTl   / Dt;
    iTl   = iTl   / Dt;
    
    // Defining cycle size:
    eCycSize = (mwSize) (eTl) + 1;
    iCycSize = (mwSize) (iTl) + 1;
    
    /* Assigning ouputs ---------------------------------------------------*/
    plhs[0] = mxCreateDoubleMatrix(simulLen, 1, mxREAL); 
    e2eI = mxGetPr(plhs[0]); 
    
    plhs[1] = mxCreateDoubleMatrix(simulLen, 1, mxREAL);
    i2eI = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateNumericMatrix(simulLen, 1, mxUINT16_CLASS, mxREAL);
    eFR = (unsigned short *) mxGetData(plhs[2]);
    
    plhs[3] = mxCreateNumericMatrix(simulLen, 1, mxUINT16_CLASS, mxREAL);
    iFR = (unsigned short *) mxGetData(plhs[3]);
    
    
    // Allocating arrays:
    tLastSP = mxCalloc(totNnrn, sizeof(double)); /* array with the last spike time of each neuron*/
    
    aX_rec  = mxCalloc(totNnrn, sizeof(double)); /* Auxiliary recurrent AMPA variables*/
    aX_ext  = mxCalloc(totNnrn, sizeof(double)); /* Auxiliary external AMPA variables*/
    gX      = mxCalloc(totNnrn, sizeof(double)); /* Auxiliary GABA variables*/
    
    aS_rec  = mxCalloc(totNnrn, sizeof(double)); /* recurrent AMPA time course*/
    aS_ext  = mxCalloc(totNnrn, sizeof(double)); /* external AMPA time course*/
    gS      = mxCalloc(totNnrn, sizeof(double)); /* GABA time course*/
    
    V       = mxCalloc(totNnrn, sizeof(double)); /* membrane potential*/
    Ncon    = mxCalloc(totNnrn, sizeof(mwSize)); /* array with the number of outgoing connections of each neuron*/
    
    // Allocating matrices:
    /* array with the number of recurrent excitatory spikes arriving on each neuron
     * in each step of the cycle*/
    eSPC    = mxCalloc(totNnrn * eCycSize, sizeof(double)); 
    /* array with the number of recurrent inhibitory spikes arriving on each neuron
     * in each time step of the cycle*/                                                             
    iSPC    = mxCalloc(totNnrn * iCycSize, sizeof(double)); 
    /* array with the connections of each neuron*/
    A       = mxCalloc(totNnrn * totNnrn , sizeof(mwSize)); 
    
    // Allocating matrices:
    /* array with the number of recurrent excitatory spikes arriving on each neuron
     * in each step of the cycle*/
    eSPC    = mxCalloc(totNnrn * eCycSize, sizeof(double)); 
    /* array with the number of recurrent inhibitory spikes arriving on each neuron
     * in each time step of the cycle*/                                                             
    iSPC    = mxCalloc(totNnrn * iCycSize, sizeof(double)); 
    /* array with the connections of each neuron*/
    A       = mxCalloc(totNnrn * totNnrn , sizeof(mwSize)); 
    
    
    /* Connections --------------------------------------------------------
     * To generate the random connections we make use of the seed that has
     * been provided by the user (if it has been provided), seed1, so that the same
     * network configuration can be used repeatedly.*/
        
    syntemp = mxCalloc(totNnrn, sizeof(mwSize)); /* auxiliary array to build the connections of the network*/
    RAND_MAX_double = (double) RAND_MAX;
    seed1 = -seed1; /*seed for ran1 must be a negative integer*/

    for(nrn=0; nrn<totNnrn; nrn++) { /*nrn is the postsynaptic neuron*/
        //recurrent AMPA connections entering the nrn-th neuron
        for(nrn2=0;nrn2<eNnrn;nrn2++) syntemp[nrn2]=0; /*initializing the syntemp array*/
        ip = p*eNnrn + (int)(0.5+gasdev(&seed1)*sqrt((float)(p*eNnrn))); /*number of AMPA connections entering the nrn-th neuron*/
        /*The mean number of AMPA connections entering each neuron is (p*eNnrn) */
        /*The variance of the number of AMPA connections entering each neuron is (p*eNnrn) */
        for(nrn2=0; nrn2<ip; nrn2++) {     
            kn = ran1(&seed1)*eNnrn; /*randomly selecting an excitatory presynaptic neuron*/
            if (syntemp[kn]==0) {
                syntemp[kn]=1;
                A[Ncon[kn] + kn*totNnrn] = nrn; /*connection: kn -> nrn*/
                Ncon[kn]++; /*number of outgoing connections of the kn-th neuron*/
            }
            else{   
                nrn2--;
            }
        }
        //GABA connections entering the nrn-th neuron
        for(nrn2=0;nrn2<iNnrn;nrn2++) syntemp[nrn2]=0; /*initializing the syntemp array*/
        ip = p*iNnrn + (int)(0.5+gasdev(&seed1)*sqrt((float)(p*iNnrn))); /*number of GABA connections entering the nrn-th neuron*/
        /*The mean number of GABA connections entering each neuron is (p*iNnrn) */
        /*The variance of the number of GABA connections entering each neuron is (p*iNnrn) */
        for(nrn2=0; nrn2<ip; nrn2++) {
            kn = eNnrn + ran1(&seed1)*iNnrn; /*randomly selecting an inhibitory presynaptic neuron*/     
            if (syntemp[kn - eNnrn]==0) {
                syntemp[kn - eNnrn]=1;
                A[Ncon[kn] + kn*totNnrn] = nrn; /*connection: kn -> nrn*/
                Ncon[kn]++; /*number of outgoing connections of the kn-th neuron*/
            }
            else{   
                nrn2--;
            }
        }
    }
    mxFree(syntemp);
    // check of the average number of outgoing connections ---------------------
    avNcon_outgo=0;
    sdNcon_outgo=0;
    for(nrn=0; nrn<totNnrn; nrn++) {    
        avNcon_outgo += Ncon[nrn];
    }
    avNcon_outgo = avNcon_outgo/totNnrn;
    for(nrn=0; nrn<totNnrn; nrn++) {
        sdNcon_outgo += pow(Ncon[nrn]-avNcon_outgo,2);
    }
    sdNcon_outgo = sqrt(sdNcon_outgo/(totNnrn-1));
   
    mexPrintf("Av number of outgoing connections per neuron = %i \nStd of the number of outgoing connections per neuron = %i\n", avNcon_outgo, sdNcon_outgo);
    
    //--------------------------------------------------------------------------
    
    eCounter=0;
    iCounter=0;
    
    for(nrn=0; nrn<totNnrn; nrn++) {
        // Initializing the membrane potential to V_leaky:
        V[nrn] = V_leaky;
        // Initializing time of last spike to -simulLen:
        tLastSP[nrn] = -simulLen;
    }
     
    mexPrintf("Simulation_time=%.0f ms\n", simulLen*Dt);
    
    // Initializing second random seed. This seed will affect the
    // poisson random variates.
    srand(seed2);
   
    mexPrintf("\nSimulation started.\n");    
    mexEvalString("drawnow;");
    
    // Loop over time -----------------------------------------------------
    for(t=0; t<simulLen; t++) {
        
        // Needed for poisson random variate generations 
        exp_x2eFR = exp(-x2eFR[t]);
        exp_x2iFR = exp(-x2iFR[t]);
        
        eCycIndx = t % eCycSize; // time index inside the cycle
        iCycIndx = t % iCycSize;
           
        /*Spike gets to all other neurons after a latency eTl
         * (if the pre-synaptic is excitatory) or (iTl if the
         *  pre-synaptic neuron is inhibitory).*/   
        eCycSP = (t + (mwSize) eTl) % eCycSize; // latency for excitatory spikes
        iCycSP = (t + (mwSize) iTl) % iCycSize; // latency for inhibitory spikes
        
        // Loop over neurons ----------------------------------------------
        for(nrn=0; nrn<totNnrn; nrn++) {
            
            /* Given the differential equation
             * 
             *       X'(t) = dX(t)/dt = F(X(t))
             * 
             *  we will use the following numerical method (that is the "midpoint method"):
             * 
             *       1. k1 = F(X(t)) * Dt/2
             *       2. rk = X(t + Dt/2) = X(t) + k1
             *       3. k2 = F(rk) * Dt
             *       4. X(t + Dt) = X(t) + k2
             * 
             *  Explanation: we have
             * 
             *       X(t + Dt)   = X(t) + F(X(t + Dt/2)) * Dt = X(t) + k2
             *       X(t + Dt/2) = X(t) + F(X(t)) * Dt/2
             *
             *
             *  Synaptic currents ------------------------------------------
             *
             *  In what follows, S is the synaptic current time course (s_syn in the paper) 
             *  produced by the activation of a single synapse at t=0, (see eq. 7 Cavallari et al 2014)
             * 
             *      S(t) = Tm/(Td-Tr) * [exp(-t/Td) - exp(-t/Tr)]
             * 
             *  and X is an auxiliary variable that we introduce to rewrite the equation for S(t) as 
             * 
             *       F(S(t)) = S'(t) = (X(t) - S(t)) / Td
             *       G(X(t)) = X'(t) = (C - X(t)) / Tr
             * 
             *  where C = Tm * delta(t). delta(t) is the delta function (that models the presence of a spike in t=0)
             *  Note that the function delta(t) is implemented in the code as 1/Dt
             *  To obtain the whole synaptic time course for a given cell, we have to replace
             *  delta(t) in the expression of C with the compound spike train of all presynaptic neuron connected to the cell:
             *  sum_i delta(t-t_i) 
             *
             *  Therefore, by applying the midpoint method, 
             *  we end up with the following numerical procedure:
             * 
             *       1) X_k1 = X'(t)*(Dt/2) = ((C - X(t)) / Tr) * (Dt/2)
             *          S_k1 = S'(t)*(Dt/2) = ((X(t) - S(t)) / Td) * (Dt/2)
             *
             *       2) X_rk = X(t) + X_k1
             *          S_rk = S(t) + S_k1
             * 
             *       3) X_k2 = X'(t+Dt/2)*Dt = ((C - X_rk) / Tr) * Dt
             *          S_k2 = S'(t+Dt/2)*Dt = ((X_rk - S_rk) / Td) * Dt
             * 
             *       4) X(t+1) = X(t) + X_k2
             *          S(t+1) = S(t) + S_k2
             *
             *
             *  Membrane potential -----------------------------------------
             *  
             *  The equation is:
             *     dV/dt = (V_leak - V(t) - Rm*synI(t)) / Tm , 
             *     where synI is the total synaptic current entering the neuron, I_tot
             *     
             *     We introduced for convenience the auxiliary variables:
             *     
             *     1) V_k1 = V'(t)*(Dt/2) = ((V_leak - V(t) - Rm*synI(t)) / Tm) * (Dt/2)
             *    
             *     2) V_rk = V(t) + V_k1
             *
             *     3) V_k2 = ((V_leak - V_rk - Rm*synI_rk) / Tm) * (Dt/2)
             *
             *     4) V(t+1) = V(t) + V_k2
             *
             *
             *  Note: in what follows, the Dt steps do not apear since all
             *  time constants have been normalized in Dt units.
             *
             */
            
            // Synaptic currents -----------------------------------------------
            
            
            /* AMPA synapses ----------------------------------------------*/
            if(nrn<eNnrn) {
                
             /* Poisson-distributed variates -------------------------------
             *  Generating external input spike-count according to a Poisson
             *  process with mean xFR[t]. Note that the following method for
             *  generating random numbers from a Poisson distribution is
             *  efficient only for small xFR (lower or equal to 12).
             *
             *  Knuth's algorithm for poisson distribution:
             */  
             // Number of external spikes impinging on the nrn-th excitatory neuron
                xSPC = -1;
                p = 1.0;
                do {
                    xSPC += 1.0;
                    rndNum = ((double) rand()) / RAND_MAX_double;
                    p *= rndNum;
                } while (p > exp_x2eFR);
                
                /*I_AMPA current entering the nrn-th excitatory neuron at time t*/
                aI=e2eJ*aS_rec[nrn]+x2eJ*aS_ext[nrn]; 
                
                e2eI[t] += aI; /*Sum of the I_AMPA currents on all the excitatory neurons at time t*/
                
                
                // E => E  ------------------------------------------------------
                
                // number of recurrent excitatory spikes impinging on the nrn-th excitatory neuron  
                e2eSPC = eSPC[nrn + eCycIndx*totNnrn];
                C = eTm * e2eSPC;
                
                aX_k1 = ((C - aX_rec[nrn]) / e2eTr) / 2;           
                aS_k1 = ((aX_rec[nrn] - aS_rec[nrn]) / e2eTd) / 2; 
                
                aX_rk = aX_rec[nrn] + aX_k1;  //  X(t+1/2)
                aS_rk_rec = aS_rec[nrn] + aS_k1; //  S(t+1/2)
                
                aX_k2 = ((C - aX_rk) / e2eTr);                  
                aS_k2 = ((aX_rk - aS_rk_rec) / e2eTd);    
                
                aX_rec[nrn] += aX_k2;  //  X(t+1)    
                aS_rec[nrn] += aS_k2;  //  S(t+1)    
                
         
                // Ext => E  ---------------------------------------------------------
                
                C = eTm * xSPC;
                
                aX_k1 = ((C - aX_ext[nrn]) / e2eTr) / 2;
                aS_k1 = ((aX_ext[nrn] - aS_ext[nrn]) / e2eTd) / 2;
                
                aX_rk = aX_ext[nrn] + aX_k1;  
                aS_rk_ext = aS_ext[nrn] + aS_k1;
                
                aX_k2 = ((C - aX_rk) / e2eTr);
                aS_k2 = ((aX_rk - aS_rk_ext) / e2eTd);
                
                aX_ext[nrn] += aX_k2;
                aS_ext[nrn] += aS_k2;
            }
            
            else {
                
                // Number of external spikes impinging on the nrn-th inhibitory neuron 
                xSPC = -1;
                p = 1.0;
                do {
                    xSPC += 1.0;
                    rndNum = ((double) rand()) / RAND_MAX_double;
                    p *= rndNum;
                } while (p > exp_x2iFR);
                
                /*I_AMPA current entering the nrn-th inhibitory neuron at time t*/
                aI=e2iJ*aS_rec[nrn]+x2iJ*aS_ext[nrn];
                
                
                // E => I  --------------------------------------------------------
                
                // number of recurrent excitatory spikes impinging on the nrn-th inhibitory neuron 
                e2iSPC = eSPC[nrn + eCycIndx*totNnrn]; 
                C = iTm * e2iSPC;
                
                aX_k1 = ((C - aX_rec[nrn]) / e2iTr) / 2;
                aS_k1 = ((aX_rec[nrn] - aS_rec[nrn]) / e2iTd) / 2;
                
                aX_rk = aX_rec[nrn] + aX_k1;  
                aS_rk_rec = aS_rec[nrn] + aS_k1; 
                
                aX_k2 = ((C - aX_rk) / e2iTr);
                aS_k2 = ((aX_rk - aS_rk_rec) / e2iTd);
                
                aX_rec[nrn] += aX_k2;
                aS_rec[nrn] += aS_k2;
                
                
                // Ext => I   ------------------------------------------------------
                
                C = iTm * xSPC;
                
                aX_k1 = ((C - aX_ext[nrn]) / e2iTr) / 2;
                aS_k1 = ((aX_ext[nrn] - aS_ext[nrn]) / e2iTd) / 2;
                
                aX_rk = aX_ext[nrn] + aX_k1;
                aS_rk_ext = aS_ext[nrn] + aS_k1;
                
                aX_k2 = ((C - aX_rk) / e2iTr);
                aS_k2 = ((aX_rk - aS_rk_ext) / e2iTd);

                aX_ext[nrn] += aX_k2;
                aS_ext[nrn] += aS_k2;
            }
            
            /* GABA synapses: ----------------------------------------------------*/
            if(nrn<eNnrn){
                
                // I => E
                
                /*I_GABA current entering the nrn-th excitatory neuron at time t*/
                gI=i2eJ*gS[nrn]; 
                
                i2eI[t] += gI; /*Sum of the I_GABA currents on all the excitatory neurons at time t*/
                
                // number of recurrent inhibitory spikes impinging on the nrn-th excitatory neuron
                i2eSPC = iSPC[nrn + iCycIndx*totNnrn];
                C = eTm * i2eSPC;
                
                gX_k1 = ((C - gX[nrn]) / iTr) / 2;
                gS_k1 = ((gX[nrn] - gS[nrn]) / iTd) / 2;
                
                gX_rk = gX[nrn] + gX_k1;
                gS_rk = gS[nrn] + gS_k1;
                
                gX_k2 = ((C - gX_rk) / iTr);
                gS_k2 = ((gX_rk - gS_rk) / iTd);
                
                gX[nrn] += gX_k2;
                gS[nrn] += gS_k2;
            }
            
            else{
                
                /*/ I => I*/
                
                /*I_GABA current entering the nrn-th inhibitory neuron at time t*/
                gI=i2iJ*gS[nrn];
                
                // number of recurrent inhibitory spikes impinging on the nrn-th inhibitory neuron
                i2iSPC = iSPC[nrn + iCycIndx*totNnrn]; 
                C = iTm * i2iSPC;
                
                gX_k1 = ((C - gX[nrn]) / iTr) / 2;
                gS_k1 = ((gX[nrn] - gS[nrn]) / iTd) / 2;
                
                gX_rk = gX[nrn] + gX_k1;
                gS_rk = gS[nrn] + gS_k1;
                
                gX_k2 = ((C - gX_rk) / iTr);
                gS_k2 = ((gX_rk - gS_rk) / iTd);
                
                gX[nrn] += gX_k2;
                gS[nrn] += gS_k2;
            }
            
            synI = aI + gI;  // Synaptic current, I_tot(t)
            
            // Membrane potential computation------------------------------------------------
            
            // Excitatory neurons
            if(nrn<eNnrn) {
                /* If the neuron is in the refractory period --------------------------------*/
                if (t-tLastSP[nrn]<eTrp) {
                    V[nrn] = eVres;
                }
                else {
                    V_k1 = (( - V[nrn] + V_leaky - eRm * synI) / eTm) / 2;
                    V_rk = V[nrn] + V_k1;
                    synI_rk = e2eJ * aS_rk_rec + x2eJ * aS_rk_ext + i2eJ * gS_rk;
                    V_k2 = (- V_rk + V_leaky - eRm * synI_rk) / eTm;  
                    V[nrn] += V_k2; /*    V(t+1)   */
                }
            }
  
            // Inhibitory neurons
            else {
                if (t-tLastSP[nrn]<iTrp) {  
                    V[nrn] = iVres;
                }
                else{
                    V_k1 = (( - V[nrn] + V_leaky - iRm * synI) / iTm) / 2;
                    V_rk = V[nrn] + V_k1; 
                    synI_rk = e2iJ * aS_rk_rec + x2iJ * aS_rk_ext + i2iJ * gS_rk;
                    V_k2 = (- V_rk + V_leaky - iRm * synI_rk) / iTm;
                    V[nrn] += V_k2;
                }
            }

            // Neuron above Vthr ------------------------------------------
            
            if(V[nrn]>Vthr) {
                tLastSP[nrn] = t;
                if(nrn<eNnrn){
                    V[nrn] = eVres;
                    eCounter++;}
                else {
                    V[nrn] = iVres;
                    iCounter++;}
                for(con=0; con<Ncon[nrn]; con++) {
                    a = A[con + nrn*totNnrn];
                    if(nrn<eNnrn)
                        eSPC[a + eCycSP * totNnrn]++;
                    else
                        iSPC[a + iCycSP * totNnrn]++;
                }
            }
            /* We have to make the current cycle free for spikes to come:*/
            eSPC[nrn + eCycIndx * totNnrn] = 0;
            iSPC[nrn + iCycIndx * totNnrn] = 0;          
        } // end of loop over neurons -----------------------------------------
        
        if (t<simulLen-1){
            eFR[t+1] = eCounter;
            eCounter=0;
            iFR[t+1] = iCounter;
            iCounter=0;
        }
    } // end of loop over time -------------------------------------------------------
    
    mxFree(tLastSP);
    mxFree(aX_rec);
    mxFree(aX_ext);
    mxFree(gX);
    mxFree(aS_rec);
    mxFree(aS_ext);
    mxFree(gS);
    mxFree(V);
    mxFree(Ncon);
    mxFree(eSPC);
    mxFree(iSPC);
    mxFree(A);

    printf("Simulation done.\n\n");  
} // main end
            
            
            
            