
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "library.h"

using namespace std;

  

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}


double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exponential(double lamda){
   if (lamda>0){ 
      double x=Rannyu();
      //SetPx(x);
      return -log(1-x)/lamda;
   }
   else{
      cerr << "lamda <0 parametro non valido" << endl;
      return 0;
   }
}
double Random :: Lorentz (double mean, double gamma){
    if (gamma>0){ 
      double x=Rannyu();
      //SetPx(x);
      return mean+gamma*tan(M_PI*(x-0.5));
   }
   else{
      cerr << "gamma <0 parametro non valido" << endl;
      return 0;
   }
}

int Random :: DisUni (int n){
   int i=1;
   double x=Rannyu(),a=0;
   while(a!=-1){
      if(a/n<=x && x<(a+1)/n){
         a=-1;
         return i;
      }
      else{
         a++;
         i++;
      }
         
   }
   return 1;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random :: MeanNVarUnif(int N){
   double mean=0;
   for(int i=0; i<N; i++){
      mean+=Rannyu()/N;
   }
   return mean;
}

double Random :: MeanNVarUnif(int N, double min, double max){
   
   double mean=0;
   for(int i=0; i<N; i++){
      mean+=Rannyu(min,max)/N;
   }
   return mean;
}

double Random :: MeanNVarExp(int N, double lamda){
   
   double mean=0;
   for(int i=0; i<N; i++){
      mean+=Exponential(lamda)/N;
   }
   return mean;
}

double Random :: MeanNVarGaus(int N, double mean, double sigma){

   double meanSum=0;
   for(int i=0; i<N; i++){
      meanSum+=Gauss(mean,sigma)/N;
   }
   return meanSum;
}
double Random :: MeanNVarLor(int N, double mean, double gamma){
   
   double meanSum=0;
   for(int i=0; i<N; i++){
      meanSum+=Lorentz(mean,gamma)/N;
   }
   return meanSum;
}


FileInt :: FileInt(){}  //Per le interazioni con i file

FileInt :: ~FileInt(){}

void FileInt :: CreateCleanFile(string namefile){ //Per creare un file vuoto o per svuotarne uno già esistente
   ofstream Stream;
   Stream.open(namefile);
   Stream.close();
   return;
} 

/*void FileInt :: WriteDouble(ofstream Stream, double x){ //scrive su un file un double
   
   if(Stream.is_open()){
      Stream << x <<  endl;
   }
   else{
      cerr << "Output stream is not open" << endl;
   }
   return;
}*/

Operation :: Operation(){}  //Per operazioni generali tipo la media

Operation :: ~Operation(){}

void Operation :: ResetMean(){ 
   mean=0;
   counter=0;
   return;
} 

void Operation :: UpdateMean(double x){
   mean+=x;
   counter+=1;
   return;
} 

double Operation :: GetMean(){
   return mean/counter;
}

void Operation :: ResetMean2(){ 
   mean2=0;
   counter2=0;
   return;
} 

void Operation :: UpdateMean2(double x){
   mean2+=x*x;
   counter2+=1;
   return;
} 

double Operation :: GetMean2(){
   return mean2/counter2;
}

double Operation ::  ErrStat(int n){
   if(counter!=counter2){
      cerr << "I contatori delle medie sono differenti" << endl;
      return 0;
   }
   if(n==0)
      return 0;
   else{
      return sqrt((GetMean2()-GetMean()*GetMean())/n);
   }
}

void Operation :: SetMinMax(double a,double b){
   if(a<=b){
      min=a;
      max=b;
   }
   else{
      min=b;
      max=a;
   }
   return;
}

void Operation :: UpdateMinMax(double x){
   if(x<min)
      min=x;
   if(x>max)
      max=x;
   return;
};
double Operation :: GetMin(){
   return min;
}
double Operation :: GetMax(){
   return max;
}
