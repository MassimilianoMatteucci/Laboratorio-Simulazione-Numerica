#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
using namespace std;

double Psi(double x, double mean,double sigma){
   return pow(M_E,-(x-mean)*(x-mean)/(2*sigma*sigma))+pow(M_E,-(x+mean)*(x+mean)/(2*sigma*sigma));
}


double V(double x){
   return pow(x,4)-5./2.*x*x;
}

double Funz(double x, double mean,double sigma){
   
   return V(x)+(-1/(sigma*sigma) +(x-mean)*(x-mean)*pow(M_E,-(x-mean)*(x-mean)/(2*sigma*sigma))/(pow(sigma,4)*Psi(x,mean,sigma))+(x+mean)*(x+mean)*pow(M_E,-(x+mean)*(x+mean)/(2*sigma*sigma))/(pow(sigma,4)*Psi(x,mean,sigma)))/-2.;
}
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
  // Stream.CreateCleanFile(filename+dat); //creo/pulisco i file di cui ho bisogno
   Op.ResetMean();
   Op.ResetMean2();
   double xstart=0.,x=0.,x1,t=1;
   double M=1000000, N=1000, n=M/N,counter=0;
   double k,a,mean=0,mu=0.80996,sigma=0.624198;
   double L=6.,dL=0.06; //Variabili per Psi^2
   double v[int(L/dL)];
   for(int i=0;i<int(L/dL);i++){v[i]=0;}
      
   // col metropolis campiono Psi^2 per ogni x calcolo Funz e ne faccio la media quello dovrebbe essere l'integrale
    

   /*for(int i=0;i<100000;i++){ //Per trovare valore di t per cui il Metropolis ha accettazione circa 50%
      if(i%100==0){
         t=0.1;
         k=0;
         counter=0;
         while((counter<49 or counter>51) and k<100){
            counter=0;
            for(int j=0;j<100;j++){
               x1=x+rnd.Rannyu(-t,t);
               Op.SetMinMax(1,Psi(x1,mu,sigma)*Psi(x1,mu,sigma)/(Psi(x,mu,sigma)*Psi(x,mu,sigma)));
               if(rnd.Rannyu()<=Op.GetMin())
                  counter++;
            }
            if (counter>52)
               t=t+0.1;
            if(counter<48)
               t=abs(t-0.1);
            k++;
         }
         //cout << k << " " << counter << " " << t<< endl;
      }
      
      x1=x+rnd.Rannyu(-t,t);
   
      Op.SetMinMax(1,Psi(x1,mu,sigma)*Psi(x1,mu,sigma)/(Psi(x,mu,sigma)*Psi(x,mu,sigma)));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         
      }
   }*/
   //cout << x <<" "<< t << " "<< k;
   
   
   double attempted=0,success=0;
// CALCOLO LINTEGRALE
  
   oStream.open(filename+dat,ios::app);
   

   Op.ResetMean();
   Op.ResetMean2();
      
   for(int i=0;i<N;i++){
      for(int j=0;j<n;j++){
      //T
         x1=x+rnd.Rannyu(-1.8,1.8);
         
      //A
         Op.SetMinMax(1,Psi(x1,mu,sigma)*Psi(x1,mu,sigma)/(Psi(x,mu,sigma)*Psi(x,mu,sigma)));
               
         if(rnd.Rannyu()<=Op.GetMin()){
            x=x1;
         }
         for(int z=0;z<int(L/dL);z++){
            if(-L/2.+dL*z<=x && x<=-L/2.+dL*(z+1)){
               v[z]++;
            }
         }
         mean+=Funz(x,mu,sigma);
      }  
      mean=mean/(double)n;
      //cout << success/attempted << endl;
         
      Op.UpdateMean(mean);
      Op.UpdateMean2(mean);
      mean=0;
      oStream<<i<<" "<< Op.GetMean()<< " " << Op.ErrStat() << endl; 
   }
        
   oStream.close();
   
   oStream.open("Psi");
   for(int i=0;i<int(L/dL);i++){
      oStream << -L/2.+dL*(i+1/2.) << " " <<v[i] <<endl;
   }
   return 0;
}