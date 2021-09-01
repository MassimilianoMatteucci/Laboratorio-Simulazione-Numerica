/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#define _USE_MATH_DEFINES
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
 
  string stato="Solid"; //per indicare in quale stato sto lavorando Solid,Liquid,Gas
  //for(int i=0;i<10;i++){
    bool restart=false;
    //if(i==0) {restart=true;};
    Input(restart);
    //SetForTemp(1.2);
    for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation{
    
      Reset(iblk);//Inizialization 
      int nconf = 1;
      for(int istep=1; istep <= nstep; ++istep){
        Move();
               //Move particles with Verlet algorithm
        if(istep*iblk==nblk*(nstep-1))
          ConfOld();
      
        if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if(istep%10 == 0){
          Measure(); 
          MeanMeasure(10);
          Accumulate();   //Properties measurement
          
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        //nconf += 1;
        }
      
      }
    G(iblk);
    }
               
    ConfFinal(); //Write final configuration to restart          
    PrintMean(stato);
  //}
  return 0;
}


void Input(bool restart){ //Prepare all stuff for the simulation

  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  if(restart==true){
//Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

//Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  else{
    cout << "Simulation starts from a previous run" << endl;
    cout << "Read initial configuration from file config.final " << endl;
    ReadConf.open("config.final");
      for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
    cout << "Read configuration at time t-dt from file old.final " << endl << endl;
    ReadConf.open("old.final");
      for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();
  }

  for(int i=0;i<3;i++){
    mean_etot[i]=0;
    mean_kin[i]=0;
    mean_temp[i]=0;
    mean_pot[i]=0;
  }
  for(int i=0;i<2;i++){
    mean2_etot[i]=0;
    mean2_kin[i]=0;
    mean2_temp[i]=0;
    mean2_pot[i]=0;
  }
   return;
}

void SetForTemp(double T){
  double t=0,Temp,f;
  Move();
  for(int i=0;i<npart; ++i){
    vx[i] = Pbc(x[i] - xold[i])/(delta);
    vy[i] = Pbc(y[i] - yold[i])/(delta);
    vz[i] = Pbc(z[i] - zold[i])/(delta);
  }
  
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  Temp= (2.0 / 3.0) * t/(double)npart;
  f=sqrt(T/Temp);
  
  for(int i=0;i<npart; ++i){
    vx[i] *= f;
    vy[i] *= f;
    vz[i] *= f;
    xold[i] = Pbc(x[i]-delta*vx[i]);
    yold[i] = Pbc(y[i]-delta*vy[i]);
    zold[i] = Pbc(z[i]-delta*vz[i]);
  }
 
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot/*, Ekin, Etot, Temp*/;

  Epot.open("output_epot.dat",ios::app);
  /*Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);*/

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
 
     bool find=false;
     // double r,step;
      //r=(double)(i-0)*box/((double)nbins*2);
     // step=box/((double)nbins*2);
      int k=0;
      while(find==false && k<nbins) {
        if((double)(k)*box/((double)nbins*2)<=dr && dr<=(double)(k)*box/((double)nbins*2)+box/((double)nbins*2)){
          walker[k]+=2.0;
          find=true;
        }
        k++;
      }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    /*Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;*/

    Epot.close();
   /* Ekin.close();
    Temp.close();
    Etot.close();*/

    return;
}

void MeanMeasure(int N){
  ofstream Epot, Ekin, Etot, Temp;
  int Nblocks;
  Nblocks=N;
  /*mean_etot[0]+=stima_etot/(nstep/(Nblocks*10));
  mean_kin[0]+=stima_kin/(nstep/(Nblocks*10));
  mean_temp[0]+=stima_temp/(nstep/(Nblocks*10));*/
  mean_pot[0]+=stima_pot/(nstep/(Nblocks*10));
  
  /*mean_etot[1]+=stima_etot/(nstep/(Nblocks*10));
  mean_kin[1]+=stima_kin/(nstep/(Nblocks*10));
  mean_temp[1]+=stima_temp/(nstep/(Nblocks*10));*/
  mean_pot[1]+=stima_pot/(nstep/(Nblocks*10));

  if(count%(nstep/(Nblocks*10))==0){
   /* mean2_etot[0]+=mean_etot[1]*mean_etot[1];
    mean2_kin[0]+=mean_kin[1]*mean_kin[1];
    mean2_temp[0]+=mean_temp[1]*mean_temp[1];*/
    mean2_pot[0]+=mean_pot[1]*mean_pot[1];
    /*mean_etot[1]=0;
    mean_kin[1]=0;
    mean_temp[1]=0;*/
    mean_pot[1]=0;
  }
  
  if(count==nstep/10){
    Epot.open("ave_epot.out",ios::app);
   /* Ekin.open("ave_ekin.out",ios::app);
    Temp.open("ave_temp.out",ios::app);
    Etot.open("ave_etot.out",ios::app);*/

    mean_pot[2]=mean_pot[0]/Nblocks;
    /*mean_kin[2]=mean_kin[0]/Nblocks;
    mean_temp[2]=mean_temp[0]/Nblocks;
    mean_etot[2]=mean_etot[0]/Nblocks;*/

    mean2_pot[1]=mean2_pot[0]/Nblocks;
    /*mean2_kin[1]=mean2_kin[0]/Nblocks;
    mean2_temp[1]=mean2_temp[0]/Nblocks;
    mean2_etot[1]=mean2_etot[0]/Nblocks;*/

    Epot << mean_pot[2] << " " << sqrt(abs(mean2_pot[1]-pow(mean_pot[2],2))/(Nblocks-1)) << endl;
    /*Ekin << mean_kin[2] << " " << sqrt(abs(mean2_kin[1]-pow(mean_kin[2],2))/(Nblocks-1))  << endl;
    Temp << mean_temp[2] << " " << sqrt(abs(mean2_temp[1]-pow(mean_temp[2],2))/(Nblocks-1)) << endl;
    Etot << mean_etot[2] << " " << sqrt(abs(mean2_etot[1]-pow(mean_etot[2],2))/(Nblocks-1)) << endl;
    */
    Epot.close();
    /*Ekin.close();
    Temp.close();
    Etot.close();*/

    /*mean_etot[0]=0;
    mean_kin[0]=0;
    mean_temp[0]=0;*/
    mean_pot[0]=0;
     
    /*mean2_etot[0]=0;
    mean2_kin[0]=0;
    mean2_temp[0]=0;*/
    mean2_pot[0]=0;
    count=0;
  }
  count++;
  return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<nbins; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
   
   
}

void G(int iblk){
  ofstream Epot, Ekin, Etot, Temp, Gofr, Gave;
  Gofr.open("output.gofr.0",ios::app);
  Gave.open("output.gave.0",ios::app);

    //g(r)
    
    
    for (int i = 0; i < nbins; i++){
      r=(double)(i)*box/((double)nbins*2);
      dr=box/((double)nbins*2);
      dV=4.0*M_PI/3.0*(pow(r+dr,3)-pow(r,3));

      glob_av[i]+=blk_av[i]/(double)blk_norm/(rho*npart*dV);
      glob_av2[i]+=pow(blk_av[i]/(double)blk_norm/(rho*npart*dV),2);
      
      Gofr  << iblk << " " << r+dr/2.0 <<" "<< blk_av[i]/(double)blk_norm/(rho*npart*dV) << endl;
      if(iblk==nblk)
        Gave  << iblk << " " << r+dr/2.0 <<" "<< glob_av[i]/(double)iblk << " "<< Error(glob_av[i],glob_av2[i],iblk) << endl;
    }

    Gofr.close();
    Gave.close();
  return;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print configuration at time t-dt to file old.final " << endl  ;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

void PrintMean(string stato){
  ofstream Epot,Ekin,Temp,Etot,Gave;

  if(stato=="Solid" or stato=="Liquid" or stato=="Gas"){
    Epot.open("Argon"+stato+"\\ave_epot"+stato+".out",ios::app);
    /*Ekin.open("Argon"+stato+"\\ave_ekin"+stato+".out",ios::app);
    Temp.open("Argon"+stato+"\\ave_temp"+stato+".out",ios::app);
    Etot.open("Argon"+stato+"\\ave_etot"+stato+".out",ios::app);*/
    Gave.open("Argon"+stato+"\\gave"+stato+".out",ios::app);

    Epot << mean_pot[2] << " " << sqrt(abs(mean2_pot[1]-pow(mean_pot[2],2))) << endl;
   /* Ekin << mean_kin[2] << " " << sqrt(abs(mean2_kin[1]-pow(mean_kin[2],2)))  << endl;
    Temp << mean_temp[2] << " " << sqrt(abs(mean2_temp[1]-pow(mean_temp[2],2))) << endl;
    Etot << mean_etot[2] << " " << sqrt(abs(mean2_etot[1]-pow(mean_etot[2],2))) << endl;*/
    
    for (int i = 0; i < nbins; i++){
      r=(double)(i)*box/((double)nbins*2);
      dr=box/((double)nbins*2);
      Gave  << nblk << " " << r+dr/2.0 <<" "<< glob_av[i]/(double)nblk << " "<< Error(glob_av[i],glob_av2[i],nblk) << endl;
    }



    Epot.close();
    /*Ekin.close();
    Temp.close();
    Etot.close();*/
    Gave.close();
  }
  return;
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
