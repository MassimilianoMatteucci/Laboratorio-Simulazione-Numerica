#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "library.h"
#include "mpi.h"
using namespace std;



int main (int argc, char *argv[]){
   
	int size,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status stat1,stat2;
	MPI_Request req;

	double tstart=MPI_Wtime();
	

   Random rnd;  //settaggio generatore numeri casuali
   int seed[4];
   int p1, p2;
   ifstream input; 
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   	if(rank==0){ input.open("seed0.in");}
	if(rank==1){input.open("seed1.in");}
	if(rank==2){input.open("seed2.in");}
	if(rank==3){input.open("seed3.in");}
	
	
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

   int nCitta=32, j=0, DimPop=1000,Nmigr=100;
   vector <Gene> G;
   vector <Cromosome> C;
   Population P(rnd);
   vector <int> Names;
   double x,y;
   vector <double> X,Y;
   ofstream S;

	
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

   Cromosome c(rnd,G);
   for(int i=0;i<DimPop;i++){
      c.NewC();
      P.AddMember(c);
   }
  
	if(rank==0){S.open("LengthsS0.dat");}
   	if(rank==3){S.open("LengthsS3.dat");}
	if(rank==2){S.open("LengthsS2.dat");}
	if(rank==1){S.open("LengthsS1.dat");}
   
   	P.CalcProb();
	P.OrdinateAll();
   	int imesg[nCitta],imesg2[nCitta];
	int itag=1,itag2=2;
	int Rk[6];
	if(rank==0){
		Rk[0] =rnd.DisUni(40)+20;
		Rk[1] =rnd.DisUni(40)+20;
		Rk[2] =rnd.DisUni(40)+20;
		Rk[3] =rnd.DisUni(40)+20;
		Rk[4] =rnd.DisUni(40)+20;
		Rk[5]=rnd.DisUni(40)+20;

	}
	MPI_Bcast(Rk,6,MPI_INTEGER,0,MPI_COMM_WORLD);
	

	
   for(int i=1;i<601;i++){
		if(i%Nmigr==0){// Scambio fra i diversi continenti
			if(rank==0){
		
				P.OrdinateAll();
				vector <Cromosome> Best=P.GetNBest(Rk[0]+Rk[1]+Rk[2]);
				Cromosome temp;
				for(int k=0;k<Best.size()-1;k++){
					int a=rnd.DisUni(Best.size()-1-k)-1;
					temp=Best[a];
					Best[a]=Best[Best.size()-1-k];
					Best[Best.size()-1-k]=temp;
				}
			
				P.ElimLastNElem(Rk[0]+Rk[1]+Rk[2]);
		
				for(int z=0;z<Rk[0];z++){
					Cromosome NewC;
					Gene NewG;
					for(int j=0;j<nCitta;j++){imesg[j]=Best[z].GetC()[j].GetName();}
					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,1,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&stat2);

					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					}
					P.AddMember(NewC);
				
				}

				for(int z=0;z<Rk[1];z++){
					Cromosome NewC;
					Gene NewG;
					for(int j=0;j<nCitta;j++){imesg[j]=Best[z+Rk[0]].GetC()[j].GetName();}
					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,2,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&stat2);
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					}
					P.AddMember(NewC);
				
				}	

				for(int z=0;z<Rk[2];z++){
					Cromosome NewC;
					Gene NewG;
					for(int j=0;j<nCitta;j++){imesg[j]=Best[z+Rk[0]+Rk[1]].GetC()[j].GetName();}
					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,3,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,3,itag2,MPI_COMM_WORLD,&stat2);

					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					}
					P.AddMember(NewC);
				
				}
			
				P.CalcProb();
				P.OrdinateAll();
				cout <<"rank "<< rank<<" ciclo n "<< i<<endl;		
				S << P.GetBest().Length() << " "<< P.Lmean() <<endl;
			}

			if(rank==1){
			
				P.OrdinateAll();
				vector <Cromosome> Best=P.GetNBest(Rk[0]+Rk[3]+Rk[4]);
				Cromosome temp;
				for(int k=0;k<Best.size()-1;k++){
					int a=rnd.DisUni(Best.size()-1-k)-1;
					temp=Best[a];
					Best[a]=Best[Best.size()-1-k];
					Best[Best.size()-1-k]=temp;
				}
				P.ElimLastNElem(Rk[0]+Rk[3]+Rk[4]);
			
				for(int z=0;z<Rk[0];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z].GetC()[j].GetName();}

					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,0,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,0,itag,MPI_COMM_WORLD,&stat1);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}
			
				for(int z=0;z<Rk[3];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg[j]=Best[z+Rk[0]].GetC()[j].GetName();}

					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,2,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&stat2);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}

				for(int z=0;z<Rk[4];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg[j]=Best[z+Rk[0]+Rk[3]].GetC()[j].GetName();}

					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,3,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,3,itag2,MPI_COMM_WORLD,&stat2);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}

				P.CalcProb();
				P.OrdinateAll();
				cout <<"rank "<< rank<<" ciclo n "<< i<<endl;
				S << P.GetBest().Length() << " "<< P.Lmean() <<endl;
			}
		
			if(rank==2){

				P.OrdinateAll();
				vector <Cromosome> Best=P.GetNBest(Rk[1]+Rk[3]+Rk[5]);
				Cromosome temp;
				for(int k=0;k<Best.size()-1;k++){
					int a=rnd.DisUni(Best.size()-1-k)-1;
					temp=Best[a];
					Best[a]=Best[Best.size()-1-k];
					Best[Best.size()-1-k]=temp;
				}
				P.ElimLastNElem(Rk[1]+Rk[3]+Rk[5]);
				for(int z=0;z<Rk[1];z++){
					Cromosome NewC;
					Gene NewG;
				
					int imesg[nCitta];
					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z].GetC()[j].GetName();}
					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,0,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,0,itag,MPI_COMM_WORLD,&stat1);
	
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				
				}

				for(int z=0;z<Rk[3];z++){
					Cromosome NewC;
					Gene NewG;
				
					int imesg[nCitta];
					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z+Rk[1]].GetC()[j].GetName();}
					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,1,itag,MPI_COMM_WORLD,&stat1);
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				
				}

				for(int z=0;z<Rk[5];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg[j]=Best[z+Rk[1]+Rk[3]].GetC()[j].GetName();}

					MPI_Isend(&imesg[0],nCitta,MPI_INTEGER,3,itag,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg2[0],nCitta,MPI_INTEGER,3,itag2,MPI_COMM_WORLD,&stat2);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg2[j]);
						NewG.SetCart(X[imesg2[j]-1],Y[imesg2[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}

				P.CalcProb();
				P.OrdinateAll();
				cout <<"rank "<< rank<<" ciclo n "<< i<<endl;			
				S << P.GetBest().Length() << " "<< P.Lmean() <<endl;
			}

			if(rank==3){
			
				P.OrdinateAll();
				vector <Cromosome> Best=P.GetNBest(Rk[2]+Rk[4]+Rk[5]);
				Cromosome temp;
				for(int k=0;k<Best.size()-1;k++){
					int a=rnd.DisUni(Best.size()-1-k)-1;
					temp=Best[a];
					Best[a]=Best[Best.size()-1-k];
					Best[Best.size()-1-k]=temp;
				}

				P.ElimLastNElem(Rk[2]+Rk[4]+Rk[5]);
				for(int z=0;z<Rk[2];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z].GetC()[j].GetName();}

					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,0,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,0,itag,MPI_COMM_WORLD,&stat1);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}	
			

				for(int z=0;z<Rk[4];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z+Rk[2]].GetC()[j].GetName();}

					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,1,itag,MPI_COMM_WORLD,&stat1);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}

				for(int z=0;z<Rk[5];z++){
					Cromosome NewC;
					Gene NewG;

					for(int j=0;j<nCitta;j++){imesg2[j]=Best[z+Rk[2]+Rk[4]].GetC()[j].GetName();}

					MPI_Isend(&imesg2[0],nCitta,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&req);
					MPI_Recv(&imesg[0],nCitta,MPI_INTEGER,2,itag,MPI_COMM_WORLD,&stat1);
				
					for(int j=0;j<nCitta;j++){
						NewG.SetName(imesg[j]);
						NewG.SetCart(X[imesg[j]-1],Y[imesg[j]-1]);
						NewC.AddGene(NewG);
					
					}
				
					P.AddMember(NewC);
				}

				P.CalcProb();
				P.OrdinateAll();
				cout <<"rank "<< rank<<" ciclo n "<< i<<endl;			
				S << P.GetBest().Length() << " "<< P.Lmean() <<endl;
			}
		
		
		
		
		}

      	P.Evolve();
    	P.CalcProb();
		P.OrdinateAll();
		S << P.GetBest().Length() << " "<< P.Lmean() <<endl;
	}

  	S.close();
	c=P.GetBest();
   
   	if(rank==0){S.open("BestPathS0.dat");}
	if(rank==1){S.open("BestPathS1.dat");}
	if(rank==2){S.open("BestPathS2.dat");}
	if(rank==3){S.open("BestPathS3.dat");}
   	for(int i=0;i<c.GetC().size();i++){
      	S << c.GetC()[i].GetName() <<" " <<c.GetC()[i].GetX() <<" "<<c.GetC()[i].GetY()<<endl; 
   	}	
   	S.close();

	
	double tend=MPI_Wtime();
	
	cout<<"rank "<<rank<<" Time " << tend-tstart << endl;
	
	MPI_Finalize();
	
   return 0;
}
