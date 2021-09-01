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

   int Nblocchi=100, Nstep=100;
   double y=0,x=0,z=0,a=1.,v[Nstep],v2[Nstep],mean[Nstep],mean2[Nstep],M=100000; //x del baricentro la y è ininfluente 
   for(int i=0;i<Nstep;i++){
      v[i]=0;
      v2[i]=0;
      mean[i]=0;
      mean2[i]=0;
   }
   
   Operation Op;
   FileInt Stream;
   ofstream oStream;
   string filename="Dati",dat=".dat";
   Stream.CreateCleanFile(filename+dat); //creo/pulisco i file di cui ho bisogno
   Stream.CreateCleanFile(filename+"2"+dat);
   double n=M/Nblocchi;
   
   for(int i=0;i<Nblocchi;i++){
      for(int w=0;w<Nstep;w++){
         v[w]=0;
         v2[w]=0;
      }
      for(int w=0;w<n;w++){
         x=0;
         y=0;
         z=0;
         for(int j=0;j<Nstep;j++){
            switch(rnd.DisUni(6)){
               case 1:
                  x-=a;
                  break;
               case 2:
                  y-=a;
                  break;
               case 3:
                  z-=a;
                  break;
               case 4:
                  x+=a;
                  break;
               case 5:
                  y+=a;
                  break;
               case 6:
                  z+=a;
                  break;
            }
            v[j]+=(x*x+y*y+z*z);
            v2[j]+=((x*x+y*y+z*z)*(x*x+y*y+z*z));
         }
      }
      
      for(int w=0;w<Nstep;w++){
         mean[w]+=v[w]/n;
         mean2[w]+=v2[w]/n;
         
      }
   }
   for(int i=0;i<Nstep;i++){
      mean[i]/=Nblocchi;
      mean2[i]/=Nblocchi;
   }
   
   oStream.open(filename+dat, ios::app);
   for(int i=0;i<Nstep;i++){
   if(oStream.is_open()){
      /* essendo dq/|q|=|n|*dx/|x|
      dove q=x^n, dq e dx l'incertezza su q ed x
      nel nostro caso calcoliamo x ovvero la media dei quadrati dei raggi
      e ne facciamo la radice della quale ci serve l'incertezza
      quindi dq=|q|*|n|*dx/|x|*/ 
         oStream << i+1 << " " << sqrt(mean[i]) << " " << 0.5*sqrt(abs(mean2[i]-mean[i]*mean[i])/(Nblocchi-1))/sqrt(abs(mean[i])) << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   } 
   oStream.close();
   
   for(int i=0;i<Nstep;i++){
      v[i]=0;
      v2[i]=0;
   }
   double Theta=0,Phi=0;
   for(int i=0;i<M;i++){
      x=0;
      y=0;
      z=0;
      for(int j=0;j<Nstep;j++){
         Theta=rnd.Rannyu(0,M_PI);
         Phi=rnd.Rannyu(0,2*M_PI);
         x+=a*cos(Phi)*sin(Theta);
         y+=a*sin(Phi)*sin(Theta);
         z+=a*cos(Theta);
         v[j]+=x*x+y*y+z*z;
         v2[j]+=(x*x+y*y+z*z)*(x*x+y*y+z*z);
      }
      
   }
   oStream.open(filename+"2"+dat, ios::app);
   for(int i=0;i<Nstep;i++){
   if(oStream.is_open()){
      /* essendo dq/|q|=|n|*dx/|x|
      dove q=x^n, dq e dx l'incertezza su q ed x
      nel nostro caso calcoliamo x ovvero la media dei quadrati dei raggi
      e ne facciamo la radice della quale ci serve l'incertezza
      quindi dq=|q|*|n|*dx/|x| */
         oStream << i << " " << sqrt(v[i]/M) << " " << 0.5*sqrt(abs(v2[i]/M-v[i]*v[i]/(M*M)))/sqrt(v[i]/M) << endl;
      }
      else{
         cerr << "Output stream is not open" << endl;
      }
   } 
   oStream.close();
   rnd.SaveSeed();
   return 0;
}


