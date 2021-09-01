#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
using namespace std;
// e tutto sbagliato per calcoli
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
   ofstream oStream,oStream2;
   string filename="Dati",dat=".dat";
   Stream.CreateCleanFile(filename+dat); //creo/pulisco i file di cui ho bisogno
   Stream.CreateCleanFile(filename+"2"+dat);
   Op.ResetMean();
   Op.ResetMean2();
   double x=1.,y=1.,z=1.,x1,y1,z1,/*A,f100,f210,T,*/t=1;
   double M=1000000, N=100, n=M/N,counter=0;
   double k,a,mean=0;
   double r=sqrt(x*x+y*y+z*z) , f=r;
   int c=1;
   oStream.open("RtoEq.dat");
   for(int i=0;i<100000;i++){ //Per raggiungere l'equilibrio
      if(i%100==0){
         t=0.1;
         k=0;
         counter=0;
         while((counter<48 or counter>52) and k<100){ //Per avere accettazione attorno al 50%
            counter=0;
            for(int j=0;j<100;j++){
               x1=x+rnd.Rannyu(-t,t);
               y1=y+rnd.Rannyu(-t,t);
               z1=z+rnd.Rannyu(-t,t);
               Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));
               if(rnd.Rannyu()<=Op.GetMin())
                  counter++;
            }
            if (counter>52)
               t=t+0.01;
            if(counter<48)
               t=abs(t-0.01);
            k++;
         }
         //cout << k << " " << counter << " " << t<< endl;
      }
      
      x1=x+rnd.Rannyu(-t,t);
      y1=y+rnd.Rannyu(-t,t);
      z1=z+rnd.Rannyu(-t,t);
   
      Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
      oStream << sqrt(x*x+y*y+z*z) << " " << t << " " << counter << endl;
   }
   oStream.close();
   k=0;
   /*while((counter<499 or counter>501) and k<1000){
      counter=0;
      for(int i=0;i<1000;i++){
            x1=x+rnd.Rannyu(-t,t);
            y1=y+rnd.Rannyu(-t,t);
            z1=z+rnd.Rannyu(-t,t);
            Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));
            if(rnd.Rannyu()<=Op.GetMin())
               counter++;
            }
      if (counter>501)
         t=t+0.1;
      if(counter<499)
         t=abs(t-0.1);
      k++;
   }
   cout << counter <<" " << t << endl;

   while(f>0.1*t and c<10000){
      c++;
      if(c%10==1){
         f=r;
      }
      x1=x+rnd.Rannyu(-t,t);
      y1=y+rnd.Rannyu(-t,t);
      z1=z+rnd.Rannyu(-t,t);
   
       Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
      if(c%10==0){
         r=sqrt(x*x+y*y+z*z);
         f=abs(f-r);
      }
   }
   cout << x << " " <<y << " " << z << " " << c << endl;*/
   //1
   oStream.open(filename+dat);
   oStream2.open("XYZ_Psi100.dat");
   for(int i=0;i<N;i++){
      for(int j=0;j<n;j++){
      //T
         x1=x+rnd.Rannyu(-t,t);
         y1=y+rnd.Rannyu(-t,t);
         z1=z+rnd.Rannyu(-t,t);
      //A
         Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));
         a=rnd.Rannyu();
         if(a<=Op.GetMin()){
            x=x1;
            y=y1;
            z=z1;
            oStream2 << x <<" "<<y<<" "<<z<<endl;
         }
         mean+=sqrt(x*x+y*y+z*z)/(double)n;
      }
      
      Op.UpdateMean(mean);
      Op.UpdateMean2(mean);
      mean=0;
      oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
      
   }

   oStream.close();
   oStream2.close();

   counter=0;
   k=0,mean=0;
   x=1.;
   y=1.;
   z=1.;
   t=0.1;
   r=sqrt(x*x+y*y+z*z);
   f=r;
   c=1;
   /*while(f>0.1*t and c<10000){
      if(c%100==0){
         t=0.1;
         k=0;
         counter=0;
         while((counter<480 or counter>520) and k<1000){
            counter=0;
            for(int i=0;i<1000;i++){
               x1=x+rnd.Rannyu(-t,t);
               y1=y+rnd.Rannyu(-t,t);
               z1=z+rnd.Rannyu(-t,t);
               Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));
               if(rnd.Rannyu()<=Op.GetMin())
                  counter++;
            }
            if (counter>520)
               t=t+0.1;
            if(counter<480)
               t=t-0.1;
            k++;
         }
         
      }
      
      x1=x+rnd.Rannyu(-t,t);
      y1=y+rnd.Rannyu(-t,t);
      z1=z+rnd.Rannyu(-t,t);
   
      Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
      if(c%100==0)
         f-=sqrt(x*x+y*y+z*z);
      c++;
   }*/

   while((counter<499 or counter>501) and k<1000){
      counter=0;
      for(int i=0;i<1000;i++){
            x1=x+rnd.Rannyu(-t,t);
            y1=y+rnd.Rannyu(-t,t);
            z1=z+rnd.Rannyu(-t,t);
            Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));
            if(rnd.Rannyu()<=Op.GetMin())
               counter++;
            }
      if (counter>501)
         t=t+0.1;
      if(counter<499)
         t=t-0.1;
      k++;
   }
   

   while(f>0.1*t and c<10000){
      c++;
      if(c%10==1){
         f=r;
      }
      x1=x+rnd.Rannyu(-t,t);
      y1=y+rnd.Rannyu(-t,t);
      z1=z+rnd.Rannyu(-t,t);
   
      Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
      if(c%10==0){
         r=sqrt(x*x+y*y+z*z);
         f=abs(f-r);
      }
   }

   
   Op.ResetMean();
   Op.ResetMean2();
   //2
   oStream.open(filename+"2"+dat);
   oStream2.open("XYZ_Psi210.dat");
   for(int i=0;i<N;i++){
      for(int j=0;j<n;j++){
      //T
         x1=x+rnd.Rannyu(-t,t);
         y1=y+rnd.Rannyu(-t,t);
         z1=z+rnd.Rannyu(-t,t);
      //A
         Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));
         
         a=rnd.Rannyu();
         if(a<=Op.GetMin()){
            x=x1;
            y=y1;
            z=z1;
            oStream2 << x << " " << y << " " << z<<endl;
         }
         mean+=sqrt(x*x+y*y+z*z)/(double)n;
         
      }
      Op.UpdateMean(mean);
      Op.UpdateMean2(mean);
      mean=0;
      oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
     
   }
   oStream.close();
   oStream2.close();

   counter=0,mean=0, x=1,y=1,z=1;
   for(int i=0;i<100000;i++){
      if(i%100==0){
         t=0.1;
         k=0;
         counter=0;
         while((counter<48 or counter>52) and k<100){
            counter=0;
            for(int j=0;j<100;j++){
               x1=x+rnd.Gauss(0,t);
               y1=y+rnd.Gauss(0,t);
               z1=z+rnd.Gauss(0,t);
               Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));
               if(rnd.Rannyu()<=Op.GetMin())
                  counter++;
            }
            if (counter>52)
               t=t+0.1;
            if(counter<48)
               t=abs(t-0.1);
            k++;
         }
         
      }
      
      x1=x+rnd.Gauss(0,t);
      y1=y+rnd.Gauss(0,t);
      z1=z+rnd.Gauss(0,t);
   
      Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
   }
   
   k=0;
   Op.ResetMean2();
   Op.ResetMean();
   oStream.open(filename+"G"+dat);
   for(int i=0;i<N;i++){
      for(int j=0;j<n;j++){
      //T
         x1=x+rnd.Gauss(0,t);
         y1=y+rnd.Gauss(0,t);
         z1=z+rnd.Gauss(0,t);
      //A
         Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)))/(pow(M_E,-sqrt(x*x+y*y+z*z))),2));
         a=rnd.Rannyu();
         if(a<=Op.GetMin()){
            x=x1;
            y=y1;
            z=z1;
         }
         mean+=sqrt(x*x+y*y+z*z)/(double)n;
      }
      
      Op.UpdateMean(mean);
      Op.UpdateMean2(mean);
      mean=0;
      oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
   }

   oStream.close();
   
   counter=0;
   k=0,mean=0;
   x=1.;
   y=1.;
   z=1.;
   t=0.1;
   r=sqrt(x*x+y*y+z*z);
   f=r;
   c=1;

   while((counter<499 or counter>501) and k<1000){
      counter=0;
      for(int i=0;i<1000;i++){
            x1=x+rnd.Gauss(0,t);
            y1=y+rnd.Gauss(0,t);
            z1=z+rnd.Gauss(0,t);
            Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));
            if(rnd.Rannyu()<=Op.GetMin())
               counter++;
            }
      if (counter>501)
         t=t+0.1;
      if(counter<499)
         t=t-0.1;
      k++;
   }
   

   while(f>0.1*t and c<10000){
      c++;
      if(c%10==1){
         f=r;
      }
      x1=x+rnd.Gauss(0,t);
      y1=y+rnd.Gauss(0,t);
      z1=z+rnd.Gauss(0,t);
   
      Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));

      a=rnd.Rannyu();
      if(a<=Op.GetMin()){
         x=x1;
         y=y1;
         z=z1;
      }
      if(c%10==0){
         r=sqrt(x*x+y*y+z*z);
         f=abs(f-r);
      }
   }
   
   Op.ResetMean();
   Op.ResetMean2();
   oStream.open(filename+"G2"+dat);
   for(int i=0;i<N;i++){
      for(int j=0;j<n;j++){
      //T
         x1=x+rnd.Gauss(0,t);
         y1=y+rnd.Gauss(0,t);
         z1=z+rnd.Gauss(0,t);
      //A
         Op.SetMinMax(1,pow((pow(M_E,-sqrt(x1*x1+y1*y1+z1*z1)/2.)*z1)/(pow(M_E,-sqrt(x*x+y*y+z*z)/2.)*z),2));
         
         a=rnd.Rannyu();
         if(a<=Op.GetMin()){
            x=x1;
            y=y1;
            z=z1;
         }
         mean+=sqrt(x*x+y*y+z*z)/(double)n;
         
      }
      Op.UpdateMean(mean);
      Op.UpdateMean2(mean);
      mean=0;
      oStream << Op.GetMean() << " " << Op.ErrStat() << endl;
   }
   oStream.close();
   rnd.SaveSeed();
   return 0;
}


