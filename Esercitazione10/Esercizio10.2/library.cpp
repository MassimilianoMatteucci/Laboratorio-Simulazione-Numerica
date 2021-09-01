
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "library.h"

using namespace std;

Population :: Population(){}

Population :: ~Population(){}

void Population :: CalcProb(){
   double tot=0;
   LC.clear();
   Prob.clear();
   for(long long unsigned int i=0;i<P.size();i++){
      LC.push_back(P[i].Length());
      tot+=1./pow(LC[i],2);
   }
   for(long long unsigned int i=0;i<P.size();i++){
      P[i].SetProb((1./pow(LC[i],2))/tot);
      Prob.push_back((1./pow(LC[i],2))/tot);
   }
}

void Population :: OrdinateAll(){
   Cromosome temp;
   double temp2,temp3;
   int min;

   for (int i=0;i<P.size()-1;i++){
      min=i;
      for(int j=i+1;j<P.size();j++){
         if(LC[j]<LC[min]){
            min=j;
         }
      }
      temp=P[min];
      P[min]=P[i];
      P[i]=temp;

      temp2=Prob[min];
      Prob[min]=Prob[i];
      Prob[i]=temp2;

      temp3=LC[min];
      LC[min]=LC[i];
      LC[i]=temp3;
   }
   Ordinated=true;
   return;
}



void Population :: Selection(){
   Pnew.push_back(P[rnd.DisIntProb(true,Prob)]);
   Pnew.push_back(P[rnd.DisIntProb(true,Prob)]);
}

void Population :: Crossover(){
   int n,j=Pnew.size()-2;
   
   long long unsigned int k=0;
   bool find;
   Cromosome Son1=Pnew[j],Son2=Pnew[j+1];
   n=rnd.DisUni(Pnew[j].GetC().size()-3);
    
   for(long long unsigned int i=1;i<Pnew[j].GetC().size();i++){
      k=n+1,find=false;
      while(k<Pnew[j].GetC().size() and find==false){
         if(Pnew[j+1].GetC()[i].GetName()==Pnew[j].GetC()[k].GetName()){
            Son1.SetGene(k,Pnew[j+1].GetC()[i]);
            find=true;
         }
         k++;
      }
      k=n+1,find=false;

      while(k<Pnew[j+1].GetC().size() and find==false){
         if(Pnew[j].GetC()[i].GetName()==Pnew[j+1].GetC()[k].GetName()){
            Son2.SetGene(k,Pnew[j].GetC()[i]);
            find=true;
         }
         k++;
      }
   }

   if(Son1.Length()<Pnew[j].Length())
      Pnew[j]=Son1;
   if(Son2.Length()<Pnew[j+1].Length())
      Pnew[j+1]=Son2;
   /*k=0;
   for(long long unsigned int i=n+1;i<Pnew[j].GetC().size();i++){
      Pnew[j+2].SetGene(i,Son1[k]);
      Pnew[j+3].SetGene(i,Son2[k]);
      k++;
   } */
   /*swap_ranges(Pnew[j+1].GetC().begin()+n+1,Pnew[j+1].GetC().end(),G1.begin());
   swap_ranges(Pnew[j].GetC().begin()+n+1,Pnew[j].GetC().end(),G.begin());*/
   
}

void Population :: Mutation(){
   int j=Pnew.size()-2;
   Cromosome C=Pnew[j];
   Gene G;
   if(rnd.Rannyu()<0.1){  
      int n=rnd.DisUni(C.GetC().size()-2);
      G=C.GetC()[n];
      C.SetGene(n,C.GetC()[n+1]); 
      C.SetGene(n+1,G);
      Pnew.push_back(C);
   }
   C=Pnew[j];
   if(rnd.Rannyu()<0.1){  
      int m=rnd.DisUni(C.GetC().size()-1);
      int n=rnd.DisUni(C.GetC().size()-3);
      
      for(int i=0;i<n;i++){
         for(int j=m;j>0;j--){
            G=C.GetC()[Pbc(i+j)];
            C.SetGene(Pbc(i+j),C.GetC()[Pbc(i+j+1)]);
            C.SetGene(Pbc(i+j+1),G);
         }
      }
      Pnew.push_back(C);
   }
   C=P[j];
   if(rnd.Rannyu()<0.1){  
      int m=rnd.DisUni(C.GetC().size()/2-1);
      int start=rnd.DisUni(C.GetC().size()-2*m);
      
      for(int j=0;j<m;j++){
         G=C.GetC()[Pbc(start+j)];
         C.SetGene(Pbc(start+j),C.GetC()[Pbc(start+j+m)]);
         C.SetGene(Pbc(start+j+m),G);
      }

      
      Pnew.push_back(C);
   }
   C=P[j];
   if(rnd.Rannyu()<0.1){  
      int m=rnd.DisUni(C.GetC().size()-1);
      int start=rnd.DisUni(C.GetC().size()-m);
      
      for(int j=0;j<m/2;j++){
         G=C.GetC()[Pbc(start+j)];
         C.SetGene(Pbc(start+j),C.GetC()[Pbc(start-j+m-1)]);
         C.SetGene(Pbc(start-j+m-1),G);
      }
   
      Pnew.push_back(C);
   }
}

void Population :: Evolve(){
   bool mutate=false;
   while(Pnew.size()<P.size()){
      Selection();
   
      if(rnd.Rannyu()<0.9){
         Crossover();
      }
      Mutation();
         
   }
   if(Pnew.size()>P.size()){
      Pnew.resize(P.size());
      
   }
   P.swap(Pnew);
   
   Pnew.clear();
}

double Population :: Lmean(){
   double mean=0;
   int Dim=P.size();
   for(int i=0;i<Dim/2;i++){
      mean+=LC[i];
   }
   return mean/(double)(Dim/2);
}

int Population :: Pbc(int i){
   if(i>=P[0].GetC().size())
      return i-P[0].GetC().size()+1;
   return i;
}

Cromosome Population :: GetBest(){
   Cromosome Best=P[0];
   for(long long unsigned int i=1;i<P.size();i++){
      if (Best.Length()>P[i].Length())
         Best=P[i];
   }
   return Best;
}

vector<Cromosome> Population :: GetNBest(int n){
   
   int i=0;
   vector<Cromosome> Best;
   if(Ordinated==false){
      cout << "Bisogna utilizzare prima il metodo Ordinate sulla popolazione per utilizzare correttamente GetNBest" << endl;
      
   }
   Best.push_back(P[i]);
   if(n==1)
      return Best;
   while(i<P.size() && Best.size()<n){
      if (P[i].GetProb()!=P[i+1].GetProb()){
         //cout << P[i].GetProb()-P[i+1].GetProb()<<endl;
         Best.push_back(P[i+1]);
      }
      i++;
   }
   if(Best.size()<n)
      cout << "Non ci sono abbastanza cromosomi differenti per restituirne il numero desiderato"<<endl;
   return Best;
}



Cromosome :: Cromosome(){}
Cromosome :: Cromosome(Random r, vector<Gene> G ){
   rnd=r;
   Cstart=G;
}
Cromosome :: ~Cromosome(){}

double Cromosome :: Length(){
   double L=0;
   for(int i=0;i<C.size();i++){
      L+=pow(C[i].GetX()-C[Pbc(i+1)].GetX(),2)+pow(C[i].GetY()-C[Pbc(i+1)].GetY(),2);
   }
   return L;
}

int Cromosome :: Pbc(int i){
   if(i>=C.size())
      return 0;
   return i;
}

void Cromosome :: NewC(){
   vector<int> Perm;
   double P;
   C.clear();
   C.push_back(Cstart[0]);
   for(int i=1;i<Cstart.size();i++){
      Perm.push_back(i);
   }
   while(C.size()!=Cstart.size()){
      P=rnd.DisUni(Perm.size())-1;
      C.push_back(Cstart[Perm[P]]);
      Perm[P]=Perm[Perm.size()-1];
      Perm.pop_back();
   }
}

Gene :: Gene(){}

Gene :: ~Gene(){}

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

int Random :: DisIntProb (bool norm,vector<double> Prb){
   bool find=false;
   int i=0;
   double x=Rannyu(),tot=0,PrbSum=0;
   if(norm==false){
      for(int j=0;j<Prb.size();j++){
         tot+=Prb[j];
      }
   }
   else{
      tot=1.;
   }
   //cout << x << endl;
   while(find==false && i<Prb.size()){
   
      if(PrbSum/tot<=x && x<(Prb[i]+PrbSum)/tot){
         find=true;
         return i;
      }
      else{
         PrbSum+=Prb[i];
         i++;
         
      }
         
   }
   //cerr << "Non è stato estratto nessun valore utile" << endl;
   return -1;
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
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
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

double Operation ::  ErrStat(){
   if(counter!=counter2){
      cerr << "I contatori delle medie sono differenti" << endl;
      return 0;
   }
   if(counter==1)
      return 0;
   else{
      return sqrt(abs(GetMean2()-GetMean()*GetMean())/(counter-1));
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

 


