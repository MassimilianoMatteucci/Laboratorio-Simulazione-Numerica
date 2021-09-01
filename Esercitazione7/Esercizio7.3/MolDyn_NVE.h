/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
using namespace std;
#include <string> 

//parameters, observables
const int m_props=4, nbins=100;
int n_props,count=1;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;
double walker[nbins];
double dr ,dV,r;

// averages
double acc,att;
double blk_av[nbins],blk_norm,accepted,attempted;
double glob_av[nbins],glob_av2[nbins];
double mean_pot[3], mean_kin[3], mean_etot[3], mean_temp[3];
double mean2_pot[2], mean2_kin[2], mean2_etot[2], mean2_temp[2];
//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
int nblk;

//functions
void Input(bool);
void Reset(int);
void Accumulate(void);
void Move(void);
void SetForTemp(double);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
void MeanMeasure(int);
void G(int);
double Force(int, int);
void PrintMean(string);
double Pbc(double);
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
