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
   double M=100000, N=10, n=M/N,counter=0;
   double k,a,mean=0,mu=1,sigma=1;
  
   // col metropolis campiono Psi^2 per ogni x calcolo Funz e ne faccio la media quello dovrebbe essere l'integrale
    

   /*for(int i=0;i<100000;i++){
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
   xstart=x;
   double Hmin;
   double attempted=0,success=0;
// CALCOLO INTEGRALE
   double b=1,sigma1=sigma,mu1=mu,H,Hnew;
   ofstream oStream2,oStream3;
   oStream.open(filename+dat,ios::app);
   oStream2.open(filename+"2"+dat,ios::app);
   
   for(double T=10;T>0.1;T-=0.1){
      attempted=0;
      success=0;
      k=0;
      b=pow(M_E,-1./T);
     
      while(k<101){
         if(k!=0){
            do{
               /*if(sigma-b>0)
                  sigma1=sigma+rnd.Rannyu(-b,b);
               else{*/
                  
                  //sigma1=sigma+rnd.Rannyu(-sigma+0.01,sigma);
               //}
               sigma1=sigma+rnd.Rannyu(-b,b);
               mu1=mu+rnd.Rannyu(-b,b);
               counter+=1;
               if(counter==100){
                  cerr << "I valori estratti di sigma e mu non sono accettabili cambiarne il punto di partenza o il passo"<< endl;
                  return -1;
               }

            }while(Psi(x,mu1,sigma1)==0 );
            counter=0;
         }

         Op.ResetMean();
         Op.ResetMean2();
      
         for(int i=0;i<N;i++){
            for(int j=0;j<n;j++){
      //T
               x1=x+rnd.Rannyu(-1.8,1.8);
               Op.SetMinMax(1,Psi(x1,mu1,sigma1)*Psi(x1,mu1,sigma1)/(Psi(x,mu1,sigma1)*Psi(x,mu1,sigma1)));
               if(rnd.Rannyu()<=Op.GetMin()){
                  x=x1;
               }
               mean+=Funz(x,mu1,sigma1);
            }
            mean=mean/(double)n;
         //cout << success/attempted << endl;
         
            Op.UpdateMean(mean);
            Op.UpdateMean2(mean);
            mean=0;
            
         }
         if(k==0 && T==10.){
            cout << "Primo" << endl;
            H=Op.GetMean();
            Hmin=H;
            oStream2 << mu << " " << sigma << " " << Hmin<<" "<<" "<<Op.ErrStat()<<endl;
            oStream << mu << " " << sigma << " " << H<<" "<<Op.ErrStat()<<endl;
         }
         else{
            Hnew=Op.GetMean();
            if(Hnew>H){
               Op.SetMinMax(1,pow(M_E,-(Hnew-H)/double(T)));
               a=rnd.Rannyu();
               
               if(a<=Op.GetMin()){
                  mu=mu1;
                  sigma=sigma1;
                  H=Hnew;
                  oStream << mu << " " << sigma << " " << H<<" "<<Op.ErrStat()<<endl;
                  success+=1;
               }
            
            }
            else{
               mu=mu1;
               sigma=sigma1;
               H=Hnew;
               oStream << mu << " " << sigma << " " << H<<" "<<Op.ErrStat()<<endl;
               
               success+=1;
            }
            attempted+=1;

            if(Hnew<=Hmin){
               oStream2 << mu << " " << sigma << " " << Hmin<<" "<<Hnew<<" "<<Op.ErrStat()<<endl;
               Hmin=Hnew;
               
            }
         }
         k++;
      
      
      }
      cout << T<<" "<< double(success/attempted) <<" "<< success<<" "<<attempted<< endl;
   }
   oStream.close();
   oStream2.close();
   
   return 0;
}