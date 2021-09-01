#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
#include "time.h"
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

   int nCitta=32, j=0, DimPop=1000;
   vector <Gene> G,G2;
   Population P(rnd);
   Population P2(rnd);
   vector <int> Names;
   double r=1.,theta,x,y;
   vector <double> R,Theta,X,Y;
   ofstream S;
   /*S.open("CityOnCirc.dat"); //Creare le città sulla circonferenza
   for(int i=0;i<nCitta;i++){
      S << i+1 <<" "<< 1. << " " << rnd.Rannyu(0,2.*M_PI) << endl;
   }
   S.close();*/

   /*S.open("CityInSquare.dat");  //Creare le città sel quadrato
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
  
  /* for(int i=0;i<DimPop;i++){
      for(int j=0;j<Names.size();j++){
         cout << P.GetMemeber(i).GetC()[j].GetName() <<" " ;
      }
      cout << endl;
   }*/
   S.open("Lengths.dat");
   P.CalcProb();
   
   for(int i=0;i<600;i++){
      P.Evolve();
      cout << i << endl;
      P.CalcProb();
      P.OrdinateAll();
      S << P.GetBest().Length() << " "<< P.Lmean() <<endl;

   }
   S.close();

   c=P.GetBest();
   
   
   S.open("BestPath.dat");
   for(int i=0;i<c.GetC().size();i++){
      S << c.GetC()[i].GetName() <<" " <<c.GetC()[i].GetX() <<" "<<c.GetC()[i].GetY()<<endl; 
   }
   S.close();

   clock_t start,end;
   start=clock();
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
   
   G2.resize(Names.size());
   
   for(int i=0;i<Names.size();i++){
      G2[i].SetName(Names[i]);
      G2[i].SetCart(X[i],Y[i]);
   }

   Cromosome c2(rnd,G2);
   for(int i=0;i<DimPop;i++){
      c2.NewC();
      P2.AddMember(c2);
   }
  
   S.open("LengthsS.dat");
   P2.CalcProb();
   
   for(int i=0;i<800;i++){
      P2.Evolve();
      cout << i << endl;
      P2.CalcProb();
      P2.OrdinateAll();
      S << P2.GetBest().Length() << " "<< P2.Lmean() <<endl;
 
   }
   S.close();

   c2=P2.GetBest();
   
   
   S.open("BestPathS.dat");
   for(int i=0;i<c2.GetC().size();i++){
      S << c2.GetC()[i].GetName() <<" " <<c2.GetC()[i].GetX() <<" "<<c2.GetC()[i].GetY()<<endl; 
   }
   S.close();
   end=clock();

   cout <<"Tempo esecuzione GA per le citta nel quadrato " <<((double)(end-start))/CLOCKS_PER_SEC;
   return 0;
}