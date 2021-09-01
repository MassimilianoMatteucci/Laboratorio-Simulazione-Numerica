#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
using namespace std;
 
int main (int argc, char *argv[]){
   if (argc!=2){
      cerr << "Use of " << argv[0] << " <int N>" << endl;
      return -1;
   }
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

   int N=atoi(argv[1]);
   int  n=10000;
   FileInt Stream;
   ofstream oStream;
   string filename, a="DatiN=", b=to_string(N) ,c=".dat";
   filename=a+b+c;
   Stream.CreateCleanFile(filename); //creo/pulisco i file di cui ho bisogno
   oStream.open(filename, ios::app);
      for(int i=0; i<n; i++){
         if(oStream.is_open()){
            oStream << rnd.MeanNVarUnif(N, -1., 1.) << " " << rnd.MeanNVarExp(N, 1.) << " " << rnd.MeanNVarLor(N, 0, 1.)  << endl;
         }
         else{
            cerr << "Output stream is not open" << endl;
         } 
      }
   oStream.close();
   
   
   rnd.SaveSeed();
   return 0;
}


