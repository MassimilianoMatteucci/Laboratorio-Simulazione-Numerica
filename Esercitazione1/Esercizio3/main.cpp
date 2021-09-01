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

   int M=1000000, N=100, HitCount=0;
   double xb,y,x; //x del baricentro la y è ininfluente
   double l=0.5,d=1,cosTheta; // lunghezza spilli e distanza linee
   //Posso considerare solo due linee di più sarebbero ridondanti
   if(l>d)
      cerr << "La distanza fra le linee deve essere maggiore della lunghezza degli spilli" << endl;
   double n=M/N;
   Operation Op;
   FileInt Stream;
   ofstream oStream;
   string filename="Dati.dat";
   Stream.CreateCleanFile(filename); //creo/pulisco i file di cui ho bisogno
   Op.ResetMean();
   Op.ResetMean2();
   oStream.open(filename, ios::app);
   for(int i=0;i<N;i++){
      HitCount=0;
      for(int j=0;j<n;j++){
         xb=rnd.Rannyu(0,d);
         y=rnd.Rannyu();
         x=rnd.Rannyu();
         if(x*x+y*y<1){
            cosTheta=(x/((sqrt(x*x+y*y))));
            if(xb>d-xb)
               xb=d-xb;
            if(xb <=abs(l*0.5*cosTheta /*cos(Theta))*/))
               HitCount++;
         }
         else{
            j--;
         }
      }
      
      Op.UpdateMean(2*l*n/(HitCount*d));
      Op.UpdateMean2(2*l*n/(HitCount*d));
      if(oStream.is_open()){
         oStream << Op.GetMean() << " " << Op.ErrStat(i) << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      } 
   }
   oStream.close();
   rnd.SaveSeed();
   return 0;
}


