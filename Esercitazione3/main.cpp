#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
using namespace std;
// devo calcolare le C 
int main (int argc, char *argv[]){
   
   Random rnd;  //settaggio generatore numeri casuali
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl; // fine settaggio
   
   Operation Op;
   FileInt Stream;
   ofstream oStream;
   string filename="Dati",dat=".dat";
   Stream.CreateCleanFile(filename+dat); //creo/pulisco i file di cui ho bisogno
   Stream.CreateCleanFile(filename+"2"+dat);
   Op.ResetMean();
   Op.ResetMean2();
   double S0=100, K=100, r=0.1, T=1, t=0, sigma=0.25;
   double M=100000, N=100, n=M/N;
   double S=0,C=0,P=0;
   oStream.open(filename+"C"+dat, ios::app);
   for(int i=0;i<N;i++){
      C=0;
      for(int j=0;j<n;j++){
         S=S0*pow(M_E,((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,T)));
         
         Op.SetMinMax(0,S-K);
         C+=pow(M_E,-r*T)*Op.GetMax();
      }
      Op.UpdateMean(C/n);
      Op.UpdateMean2(C/n);
      if(oStream.is_open()){
         oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   }
   oStream.close();
   Op.ResetMean();
   Op.ResetMean2();

   oStream.open(filename+"P"+dat, ios::app);
   for(int i=0;i<N;i++){
      P=0;
      for(int j=0;j<n;j++){
         S=S0*pow(M_E,((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,T)));
         
         Op.SetMinMax(0,K-S);
         P+=pow(M_E,-r*T)*Op.GetMax();
      }
      Op.UpdateMean(P/n);
      Op.UpdateMean2(P/n);
      if(oStream.is_open()){
         oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   }
   oStream.close();
   Op.ResetMean();
   Op.ResetMean2();
   oStream.open(filename+"C2"+dat, ios::app);
   for(int i=0;i<N;i++){
      C=0;
      for(int j=0;j<n;j++){
         S=S0;
         for(int z=0;z<100;z++){
            S*=pow(M_E,((r-0.5*sigma*sigma)*(T/100.)+sigma*rnd.Gauss(0,1)*pow(T/100.,0.5)));
         }
         Op.SetMinMax(0,S-K);
         C+=pow(M_E,-r*T)*Op.GetMax();
         Op.SetMinMax(0,K-S);
         P+=pow(M_E,-r*T)*Op.GetMax();
      }
      Op.UpdateMean(C/n);
      Op.UpdateMean2(C/n);
      if(oStream.is_open()){
         oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   }
   oStream.close();

   Op.ResetMean();
   Op.ResetMean2();
   oStream.open(filename+"P2"+dat, ios::app);
   for(int i=0;i<N;i++){
      P=0;
      for(int j=0;j<n;j++){
         S=S0;
         for(int z=0;z<100;z++){
            S*=pow(M_E,((r-0.5*sigma*sigma)*(T/100.)+sigma*rnd.Gauss(0,1)*pow(T/100.,0.5)));
         }
         Op.SetMinMax(0,S-K);
         C+=pow(M_E,-r*T)*Op.GetMax();
         Op.SetMinMax(0,K-S);
         P+=pow(M_E,-r*T)*Op.GetMax();
      }
      Op.UpdateMean(P/n);
      Op.UpdateMean2(P/n);
      if(oStream.is_open()){
         oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   }
   oStream.close();
   rnd.SaveSeed();
   return 0;
}


