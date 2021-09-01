#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   Error Er;
   int M=100000, N=100;
   double Err=0, sum=0, Val=0, media=0, media2=0, sumVal=0, sumVal2=0 , L=M/N, y=0;
   ofstream Dati;
   Dati.open("Dati.dat");
   for(int i=0; i<N; i++){
      sum=0;
      for(int j=0; j<L; j++){
         y=rnd.Rannyu();
         sum+= y;
      }
      Val=sum/L;
      sumVal+=Val;
      sumVal2+=Val*Val;
   
      media=sumVal/(i+1);
      media2=sumVal2/(i+1);
      
      Err=Er.ErrStat(media, media2, i);
      if(Dati.is_open()){
         Dati << media << " " << Err << endl;
      }
      else{
         cerr << "Error impossible to open ofstream" << endl;
      }
      
   }
   Dati.close();
   rnd.SaveSeed();
   return 0;
}


