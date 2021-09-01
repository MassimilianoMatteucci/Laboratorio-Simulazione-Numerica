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

   int nCitta=32, j=0, DimPop=1;
   vector <Gene> G;
   vector <Cromosome> C,C2;
   Population P(rnd),P2(rnd);
   vector <int> Names;
   double r,theta,x,y;
   vector <double> R,Theta,X,Y;
   ofstream S;
   /*S.open("CityOnCirc.dat");
   for(int i=0;i<nCitta;i++){
      S << i+1 <<" "<< 1. << " " << rnd.Rannyu(0,2.*M_PI) << endl;
   }
   S.close();

  /* S.open("CityInSquare.dat");
   for(int i=0;i<nCitta;i++){
      S << i+1 <<" "<< rnd.Rannyu(0.,r) << " " << rnd.Rannyu(0.,r) << endl;
   }
   S.close();*/
   
  /* S.open("CityOnCirc.dat");
   for(int i=0;i<nCitta;i++){
      S << i+1 <<" "<< 1. << " " << 2.*M_PI/32.*i << endl;
   }
   S.close();*/

   input.open("CityOnCirc.dat");
   if (input.is_open()){
      while ( !input.eof() ){
         input >> j >> r >> theta;
         Names.push_back(j);
         R.push_back(r);
         Theta.push_back(theta);
      }
      input.close();
   }
   
   if(R[R.size()-1]==R[R.size()-2] and Theta[Theta.size()-1]==Theta[Theta.size()-2]){
      R.pop_back();
      Theta.pop_back();
      Names.pop_back();
   }
   
   G.resize(Names.size());
   
   for(int i=0;i<Names.size();i++){
      G[i].SetName(Names[i]);
      G[i].SetPol(Theta[i],R[i]);
      G[i].CarfromPol();
   }

   Cromosome c(rnd,G);
   for(int i=0;i<DimPop;i++){
      c.NewC();
      P.AddMember(c);
   }
  int k=0;
   double T=20.,beta;
   S.open("Lengths.dat");
   
   while(T>0.1){
      T-=0.1;
      for(int j=0;j<100;j++){
         P.Mutation1(0);
         beta=1/T;
         if(P.GetMemberN(0).Length()<P.GetMember(0).Length()){
            P.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P.GetMemberN(0).Length(),P.GetMember(0)).Length())){
               P.SwapNewOld();
            }
         }
         P.ClearNew();
         P.Mutation2(0);
         if(P.GetMemberN(0).Length()<P.GetMember(0).Length()){
            P.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P.GetMemberN(0).Length(),P.GetMember(0)).Length())){
               P.SwapNewOld();
            }
         }
         P.ClearNew();
         P.Mutation3(0);
         
         if(P.GetMemberN(0).Length()<P.GetMember(0).Length()){
            P.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P.GetMemberN(0).Length(),P.GetMember(0)).Length())){
               P.SwapNewOld();
            }
         }
         P.ClearNew();
         P.Mutation4(0);
         
         if(P.GetMemberN(0).Length()<P.GetMember(0).Length()){
            P.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P.GetMemberN(0).Length(),P.GetMember(0)).Length())){
               P.SwapNewOld();
            }
         }
         P.ClearNew();
         S << P.GetMember(0).Length() << endl;
      }
   }
   S.close();

   c=P.GetBest();
   
   
   S.open("BestPath.dat");
   for(int i=0;i<c.GetC().size();i++){
      S << c.GetC()[i].GetName() <<" " <<c.GetC()[i].GetX() <<" "<<c.GetC()[i].GetY()<<endl; 
   }
   S.close();

   Names.clear();
   input.open("CityInSquare.dat");
   if (input.is_open()){
      while ( !input.eof() ){
         input >> j >> x >> y;
         Names.push_back(j);
         X.push_back(x);
         Y.push_back(y);
      }
      input.close();
   }
   
   
   if(X[X.size()-1]==X[X.size()-2] and Y[Y.size()-1]==Y[Y.size()-2]){
      X.pop_back();
      Y.pop_back();
      Names.pop_back();
   }
   
   G.resize(Names.size());
   
   for(int i=0;i<Names.size();i++){
      G[i].SetName(Names[i]);
      G[i].SetCart(X[i],Y[i]);
   }

   Cromosome c2(rnd,G);
   for(int i=0;i<DimPop;i++){
      c2.NewC();
      P2.AddMember(c2);
   }
  
   S.open("LengthsS.dat");
   T=2.0,beta;
    
   while(T>0.01){
      T-=0.01;
      cout << "Temperatura " << T <<endl;
      for(int j=0;j<1200;j++){
         P2.Mutation1(0);
         beta=1./T;
         
         if(P2.GetMemberN(0).Length()<P2.GetMember(0).Length()){
            P2.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P2.GetMemberN(0).Length()-P2.GetMember(0).Length()))){
               P2.SwapNewOld();
            }
         }
         P2.ClearNew();
         P2.Mutation2(0);
         if(P2.GetMemberN(0).Length()<P2.GetMember(0).Length()){
            P2.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P2.GetMemberN(0).Length()-P2.GetMember(0).Length()))){
               P2.SwapNewOld();
            }
         }
         P2.ClearNew();
         P2.Mutation3(0);
         
         if(P2.GetMemberN(0).Length()<P2.GetMember(0).Length()){
            P2.SwapNewOld();
         }
         else{
            if(rnd.Rannyu()<pow(M_E,-beta*(P2.GetMemberN(0).Length()-P2.GetMember(0).Length()))){
               P2.SwapNewOld();
            }
         }
         P2.ClearNew();
         P2.Mutation4(0);
         
         if(P2.GetMemberN(0).Length()<P2.GetMember(0).Length()){
            P2.SwapNewOld();
         }
         else{
           
            if(rnd.Rannyu()<pow(M_E,-beta*(P2.GetMemberN(0).Length()-P2.GetMember(0).Length()))){
               P2.SwapNewOld();
            }
         }
         P2.ClearNew();
         S << P2.GetMember(0).Length() << endl;
      }
   }
   S.close();
   
   c2=P2.GetBest();
   
   S.open("BestPathS.dat");
   for(int i=0;i<c2.GetC().size();i++){
      S << c2.GetC()[i].GetName() <<" " <<c2.GetC()[i].GetX() <<" "<<c2.GetC()[i].GetY()<<endl; 
   }
   S.close();
   return 0;
}