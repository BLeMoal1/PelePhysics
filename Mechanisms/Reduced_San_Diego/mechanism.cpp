#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {8,15,23,34,55,56,83,91,101,4,5,6,7,27,60,78,98,99,205,206,207,208,209,210,211,212,213,214,0,1,2,3,9,10,11,12,13,14,16,17,18,19,20,21,22,24,25,26,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,57,58,59,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,79,80,81,82,84,85,86,87,88,89,90,92,93,94,95,96,97,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,215};

// Returns 0-based map of reaction order
void GET_RMAP
(int * _rmap)
{
for (int j=0; j<NUM_REACTIONS; ++j)
{
_rmap[j] = rmap[j];
}
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int& nspec, int ki[], int nu[])
{
const int ns[NUM_GAS_REACTIONS] =
     {4,4,4,3,2,3,2,3,3,3,4,4,4,4,4,2,3,3,4,4,4,4,4,3,4,4,4,3,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,2,4,5,4,2,4,4,4,4,4,4,5,3,4,4,4,4,4,4,4,4,4,3,4,4,4,4,3,4,4,4,4,4,3,3,3,4,4,4,4,3,4,3,3,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,5,4,4,4,4,4,3,4,4,4,4,5,5,5,3,4,4,2,4,4,4,4,4,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,4,5,5,3,3,3,3,3,3,3,3,3,3,4};
const int kiv[NUM_GAS_REACTIONS*5] =
     {1,2,4,3,0,5,4,1,3,0,5,3,1,6,0,6,4,3,0,0,1,5,0,0,0,1,3,6,0,0,4,2,0,0,0,1,4,3,0,0,1,2,7,0,0,1,7,3,0,0,1,7,5,2,0,1,7,6,4,0,7,4,2,3,0,7,3,6,2,0,7,3,6,2,0,3,8,0,0,0,7,8,2,0,0,7,8,2,0,0,1,8,5,7,0,1,8,6,3,0,8,3,6,7,0,8,3,6,7,0,8,4,7,3,0,9,4,10,0,0,9,3,10,1,0,9,7,10,3,0,9,2,10,4,0,11,9,1,0,0,1,11,9,5,0,11,4,9,3,0,11,4,10,1,0,11,3,9,6,0,11,2,9,7,0,12,11,13,9,0,1,11,14,0,0,14,1,5,11,0,14,4,11,3,0,14,3,6,11,0,14,2,11,7,0,14,7,8,11,0,13,1,12,5,0,13,3,12,6,0,13,4,12,3,0,13,2,12,7,0,13,7,12,8,0,12,1,5,15,0,12,1,5,16,0,12,3,6,16,0,12,4,14,1,0,12,15,17,1,0,12,7,18,3,0,12,2,14,3,0,12,2,18,4,0,12,17,5,0,0,12,19,1,0,0,12,1,13,0,0,12,20,0,0,0,3,16,14,1,0,2,16,9,1,3,10,16,14,9,0,16,15,0,0,0,1,15,21,5,0,3,15,14,1,0,3,15,21,6,0,4,15,9,1,0,4,15,9,5,0,2,15,10,5,0,2,15,9,1,3,15,22,1,0,0,21,4,9,1,0,21,2,11,4,0,21,6,14,1,0,21,10,9,11,0,18,1,14,5,0,18,1,6,16,0,18,3,14,6,0,18,4,14,3,0,18,2,14,7,0,18,14,1,0,0,20,1,19,5,0,20,4,19,3,0,20,3,19,6,0,20,12,19,13,0,20,19,1,0,0,20,7,19,8,0,19,1,17,5,0,19,4,17,3,0,19,4,14,12,0,19,2,17,7,0,19,2,23,0,0,23,17,7,0,0,19,17,1,0,0,17,1,24,5,0,17,3,24,6,0,17,4,12,11,0,17,4,25,1,0,17,24,19,0,0,17,2,24,7,0,17,24,1,0,0,17,22,5,0,0,24,1,22,5,0,24,22,1,0,0,24,2,14,11,0,24,2,25,4,0,24,2,22,7,0,22,4,1,26,0,22,4,9,15,0,22,2,14,9,0,22,3,27,1,0,22,3,28,6,0,27,1,12,9,0,27,4,10,15,0,27,4,26,3,0,27,12,19,9,0,1,26,9,16,0,26,3,9,1,11,26,4,9,1,0,26,2,9,3,0,26,2,9,10,1,28,3,1,26,0,28,4,21,9,0,28,2,26,4,0,28,2,21,10,0,28,2,9,11,0,25,27,1,0,0,25,1,12,11,0,25,1,27,5,0,25,4,14,11,0,25,3,27,6,0,25,2,14,9,3,25,12,19,9,1,25,7,14,11,3,25,12,9,0,0,28,4,29,9,0,28,2,29,10,0,29,21,0,0,0,29,0,21,0,0,29,2,21,2,0,29,6,21,6,0,29,5,21,5,0,29,10,21,10,0,29,9,21,9,0,13,29,21,13,0,21,2,9,30,0,30,3,0,0,0,0,30,0,3,0,2,30,2,3,0,6,30,6,3,0,5,30,5,3,0,10,30,10,3,0,9,30,9,3,0,13,30,13,3,0,21,4,41,32,0,41,32,9,1,0,6,32,9,31,0,41,31,1,6,0,41,31,1,3,0,41,31,5,3,0,41,31,1,5,4,42,31,5,32,0,27,32,33,9,0,12,32,33,1,0,33,41,27,1,0,27,31,33,6,0,33,41,12,9,0,33,4,14,32,0,44,32,34,9,0,44,31,34,6,0,34,41,44,1,0,27,34,33,44,0,5,35,41,8,0,1,35,41,7,0,35,3,2,37,0,1,35,4,37,0,4,37,41,7,0,1,37,41,6,0,42,37,41,11,0,21,37,14,41,0,12,37,44,41,0,40,1,10,37,0,40,4,10,35,0,38,1,10,41,5,11,37,38,1,0,42,36,9,41,0,5,36,1,37,0,13,36,12,37,0,6,36,3,37,0,14,36,11,37,0,14,36,38,1,0,20,36,19,37,0,1,36,41,3,0,5,36,41,6,0,21,36,41,11,0,43,36,14,41,0,9,36,10,41,0,4,36,41,2,0,22,36,27,41,0,6,36,41,8,0,4,35,2,36,0,33,35,27,1,2,34,35,12,6,2,33,36,27,1,4,33,36,25,4,0,34,36,12,6,4,34,39,44,10,6,41,2,35,0,0,41,2,35,0,0,41,2,35,0,0,41,2,35,0,0,41,3,37,0,0,10,37,39,0,0,10,37,39,0,0,41,4,36,0,0,41,4,36,0,0,10,36,40,0,0,20,21,19,43,0};
const int nuv[NUM_GAS_REACTIONS*5] =
     {-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,0,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,0,0,0,-2,1,1,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,1,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-2,1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,2,1,0,-1,-1,2,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0};
if (i < 1) {
// Return max num species per reaction
nspec = 5;
} else {
if (i > NUM_GAS_REACTIONS) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*5 + j] + 1;
nu[j] = nuv[(i-1)*5 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[], amrex::Real q_f[], amrex::Real q_r[])
{
amrex::Real c[45]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 45; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, T);

// convert to chemkin units
for (int id = 0; id < 216; ++id) {
q_f[id] *= 1.0e-6;
q_r[id] *= 1.0e-6;
}
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);
// compute the Gibbs free energy
amrex::Real g_RT[45];
gibbs(g_RT, T);

amrex::Real sc_qss[1];
comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);

}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 14.007000; // N
awt[1] = 1.008000; // H
awt[2] = 15.999000; // O
awt[3] = 12.011000; // C
awt[4] = 0.000549; // E
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int kd = 5; 
// Zero ncf
for (int id = 0; id < kd * 45; ++ id) {
 ncf[id] = 0; 
}

// N2
ncf[ 0 * kd + 0 ] = 2; // N

// H
ncf[ 1 * kd + 1 ] = 1; // H

// O2
ncf[ 2 * kd + 2 ] = 2; // O

// OH
ncf[ 3 * kd + 1 ] = 1; // H
ncf[ 3 * kd + 2 ] = 1; // O

// O
ncf[ 4 * kd + 2 ] = 1; // O

// H2
ncf[ 5 * kd + 1 ] = 2; // H

// H2O
ncf[ 6 * kd + 1 ] = 2; // H
ncf[ 6 * kd + 2 ] = 1; // O

// HO2
ncf[ 7 * kd + 1 ] = 1; // H
ncf[ 7 * kd + 2 ] = 2; // O

// H2O2
ncf[ 8 * kd + 1 ] = 2; // H
ncf[ 8 * kd + 2 ] = 2; // O

// CO
ncf[ 9 * kd + 3 ] = 1; // C
ncf[ 9 * kd + 2 ] = 1; // O

// CO2
ncf[ 10 * kd + 3 ] = 1; // C
ncf[ 10 * kd + 2 ] = 2; // O

// HCO
ncf[ 11 * kd + 3 ] = 1; // C
ncf[ 11 * kd + 1 ] = 1; // H
ncf[ 11 * kd + 2 ] = 1; // O

// CH3
ncf[ 12 * kd + 3 ] = 1; // C
ncf[ 12 * kd + 1 ] = 3; // H

// CH4
ncf[ 13 * kd + 3 ] = 1; // C
ncf[ 13 * kd + 1 ] = 4; // H

// CH2O
ncf[ 14 * kd + 3 ] = 1; // C
ncf[ 14 * kd + 1 ] = 2; // H
ncf[ 14 * kd + 2 ] = 1; // O

// T-CH2
ncf[ 15 * kd + 3 ] = 1; // C
ncf[ 15 * kd + 1 ] = 2; // H

// S-CH2
ncf[ 16 * kd + 3 ] = 1; // C
ncf[ 16 * kd + 1 ] = 2; // H

// C2H4
ncf[ 17 * kd + 3 ] = 2; // C
ncf[ 17 * kd + 1 ] = 4; // H

// CH3O
ncf[ 18 * kd + 3 ] = 1; // C
ncf[ 18 * kd + 1 ] = 3; // H
ncf[ 18 * kd + 2 ] = 1; // O

// C2H5
ncf[ 19 * kd + 3 ] = 2; // C
ncf[ 19 * kd + 1 ] = 5; // H

// C2H6
ncf[ 20 * kd + 3 ] = 2; // C
ncf[ 20 * kd + 1 ] = 6; // H

// CH
ncf[ 21 * kd + 3 ] = 1; // C
ncf[ 21 * kd + 1 ] = 1; // H

// C2H2
ncf[ 22 * kd + 3 ] = 2; // C
ncf[ 22 * kd + 1 ] = 2; // H

// C2H4OOH
ncf[ 23 * kd + 3 ] = 2; // C
ncf[ 23 * kd + 1 ] = 5; // H
ncf[ 23 * kd + 2 ] = 2; // O

// C2H3
ncf[ 24 * kd + 3 ] = 2; // C
ncf[ 24 * kd + 1 ] = 3; // H

// CH2CHO
ncf[ 25 * kd + 3 ] = 2; // C
ncf[ 25 * kd + 1 ] = 3; // H
ncf[ 25 * kd + 2 ] = 1; // O

// HCCO
ncf[ 26 * kd + 3 ] = 2; // C
ncf[ 26 * kd + 1 ] = 1; // H
ncf[ 26 * kd + 2 ] = 1; // O

// CH2CO
ncf[ 27 * kd + 3 ] = 2; // C
ncf[ 27 * kd + 1 ] = 2; // H
ncf[ 27 * kd + 2 ] = 1; // O

// C2H
ncf[ 28 * kd + 3 ] = 2; // C
ncf[ 28 * kd + 1 ] = 1; // H

// CHs
ncf[ 29 * kd + 3 ] = 1; // C
ncf[ 29 * kd + 1 ] = 1; // H

// OHs
ncf[ 30 * kd + 1 ] = 1; // H
ncf[ 30 * kd + 2 ] = 1; // O

// H3Op
ncf[ 31 * kd + 4 ] = -1; // E
ncf[ 31 * kd + 1 ] = 3; // H
ncf[ 31 * kd + 2 ] = 1; // O

// HCOp
ncf[ 32 * kd + 3 ] = 1; // C
ncf[ 32 * kd + 4 ] = -1; // E
ncf[ 32 * kd + 1 ] = 1; // H
ncf[ 32 * kd + 2 ] = 1; // O

// C2H3Op
ncf[ 33 * kd + 3 ] = 2; // C
ncf[ 33 * kd + 4 ] = -1; // E
ncf[ 33 * kd + 1 ] = 3; // H
ncf[ 33 * kd + 2 ] = 1; // O

// CH5Op
ncf[ 34 * kd + 3 ] = 1; // C
ncf[ 34 * kd + 4 ] = -1; // E
ncf[ 34 * kd + 1 ] = 5; // H
ncf[ 34 * kd + 2 ] = 1; // O

// O2n
ncf[ 35 * kd + 4 ] = 1; // E
ncf[ 35 * kd + 2 ] = 2; // O

// On
ncf[ 36 * kd + 4 ] = 1; // E
ncf[ 36 * kd + 2 ] = 1; // O

// OHn
ncf[ 37 * kd + 4 ] = 1; // E
ncf[ 37 * kd + 1 ] = 1; // H
ncf[ 37 * kd + 2 ] = 1; // O

// CHO2n
ncf[ 38 * kd + 3 ] = 1; // C
ncf[ 38 * kd + 4 ] = 1; // E
ncf[ 38 * kd + 1 ] = 1; // H
ncf[ 38 * kd + 2 ] = 2; // O

// CHO3n
ncf[ 39 * kd + 3 ] = 1; // C
ncf[ 39 * kd + 4 ] = 1; // E
ncf[ 39 * kd + 1 ] = 1; // H
ncf[ 39 * kd + 2 ] = 3; // O

// CO3n
ncf[ 40 * kd + 3 ] = 1; // C
ncf[ 40 * kd + 4 ] = 1; // E
ncf[ 40 * kd + 2 ] = 3; // O

// E
ncf[ 41 * kd + 4 ] = 1; // E

// C
ncf[ 42 * kd + 3 ] = 1; // C

// CH2
ncf[ 43 * kd + 3 ] = 1; // C
ncf[ 43 * kd + 1 ] = 2; // H

// CH3OH
ncf[ 44 * kd + 3 ] = 1; // C
ncf[ 44 * kd + 1 ] = 4; // H
ncf[ 44 * kd + 2 ] = 1; // O

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(5);
ename[0] = "N";
ename[1] = "H";
ename[2] = "O";
ename[3] = "C";
ename[4] = "E";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(45);
kname[0] = "N2";
kname[1] = "H";
kname[2] = "O2";
kname[3] = "OH";
kname[4] = "O";
kname[5] = "H2";
kname[6] = "H2O";
kname[7] = "HO2";
kname[8] = "H2O2";
kname[9] = "CO";
kname[10] = "CO2";
kname[11] = "HCO";
kname[12] = "CH3";
kname[13] = "CH4";
kname[14] = "CH2O";
kname[15] = "T-CH2";
kname[16] = "S-CH2";
kname[17] = "C2H4";
kname[18] = "CH3O";
kname[19] = "C2H5";
kname[20] = "C2H6";
kname[21] = "CH";
kname[22] = "C2H2";
kname[23] = "C2H4OOH";
kname[24] = "C2H3";
kname[25] = "CH2CHO";
kname[26] = "HCCO";
kname[27] = "CH2CO";
kname[28] = "C2H";
kname[29] = "CHs";
kname[30] = "OHs";
kname[31] = "H3Op";
kname[32] = "HCOp";
kname[33] = "C2H3Op";
kname[34] = "CH5Op";
kname[35] = "O2n";
kname[36] = "On";
kname[37] = "OHn";
kname[38] = "CHO2n";
kname[39] = "CHO3n";
kname[40] = "CO3n";
kname[41] = "E";
kname[42] = "C";
kname[43] = "CH2";
kname[44] = "CH3OH";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<46; k++) {
for (int l=0; l<46; l++) {
if(Jac[ 46 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<46; k++) {
for (int l=0; l<46; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 46 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<46; k++) {
for (int l=0; l<46; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 46 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 46;
int offset_col = nc * 46;
for (int k=0; k<46; k++) {
for (int l=0; l<46; l++) {
if(Jac[46*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 46;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 46;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 46;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 46;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<46; k++) {
for (int l=0; l<46; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 46*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[46*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 46*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,2116> Jac = {0.0};
amrex::GpuArray<amrex::Real,45> conc = {0.0};
for (int n=0; n<45; n++) {
    conc[n] = 1.0/ 45.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<46; l++) {
for (int k=0; k<46; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[46*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
