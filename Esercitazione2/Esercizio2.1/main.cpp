#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
using namespace std;
 
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

   int M=100000, N=200;
   double y,x; //x del baricentro la y è ininfluente 
   
   double n=M/N;
   Operation Op;
   FileInt Stream;
   ofstream oStream;
   string filename="Dati",dat=".dat";
   Stream.CreateCleanFile(filename+dat); //creo/pulisco i file di cui ho bisogno
   Stream.CreateCleanFile(filename+"2"+dat);
   Op.ResetMean();
   Op.ResetMean2();

   oStream.open(filename+dat, ios::app);
   for(int i=0;i<N;i++){
      y=0;
      for(int j=0;j<n;j++){
         y+=M_PI/2.*cos(M_PI/2.*rnd.Rannyu());
         
      }
      Op.UpdateMean(y/n);
      Op.UpdateMean2(y/n);
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
   oStream.open(filename+"2"+dat, ios::app);
   double z=0;
   for(int i=0;i<N;i++){
      y=0;
      for(int j=0;j<n;j++){
         x=rnd.Rannyu();
         z=rnd.Rannyu();
         if(z<=(1.-x*x/2)){
            y+=5./6.*M_PI/2*cos(M_PI/2*x)/(1.-x*x/2);
            
         }
         else{
            j--;
         }
      }
      Op.UpdateMean(y/n);
      Op.UpdateMean2(y/n);
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


