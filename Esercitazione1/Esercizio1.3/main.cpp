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

   int M=100, N=10000;
   double n[M], chi2=0, y=0;
   ofstream Dati("Dati.dat");
   
   for(int i=0; i<100; i++){
      for(int q=0; q<M; q++){
         n[q]=0;
      }
      chi2=0;
      for(int j=0; j<N; j++){ 
         y=rnd.Rannyu();
         for(int z=0; z<M; z++){ //ciclo per 'posizionare' il numero nell'intervallo corretto
            if(y>(double) z/M && y<(double) z/M +1./M)
               n[z]+=1; // ho un punto in piu nell'i-esimo intervallo
         }
      }
      for(int q=0; q<M; q++){
         chi2+= (n[q]-(double) N/M)*(n[q]-(double) N/M)*((double)M/N);
      }
      Dati.open("Dati.dat", ios :: app);
      if(Dati.is_open()){
         Dati << chi2 << endl;
      }
      else{
         cerr << "Error impossible open ofstream" << endl;
      }
      Dati.close();
   }
   
   rnd.SaveSeed();
   return 0;
}


