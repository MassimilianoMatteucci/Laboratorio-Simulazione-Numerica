

#ifndef __Library__
#define __Library__
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
using namespace std;
 
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max); //genera variabili con distribuzione uniforma
  double Gauss(double mean, double sigma); 
  double Exponential(double lamda);
  double Lorentz (double mean, double gamma);
  double MeanNVarUnif(int N); //media tra N variabili generate con distribuzione uniforme [0,1)
  double MeanNVarUnif(int N, double min, double max); //media tra N variabili generate con distribuzione uniforme [a,b)
  double MeanNVarGaus(int N, double mean, double sigma);
  double MeanNVarExp(int N, double lamda);
  double MeanNVarLor(int N, double mean, double gamma);
  int DisUni (int n);
  int DisIntProb (bool,vector<double> v);
};

class Gene{
private:
  int Name;
  double Cart[2],Pol[2];
public:
  Gene();
  ~Gene();

  void SetName(int i){Name=i;}
  void SetX(double x){Cart[0]=x;}
  void SetY(double y){Cart[1]=y;}
  void SetTheta(double theta){Pol[0]=theta;}
  void SetR(double r){Pol[1]=r;}
  void SetCart(double x,double y){Cart[0]=x; Cart[1]=y;}
  void SetPol(double theta,double r){Pol[0]=theta; Pol[1]=r;}

  double GetX(){return Cart[0];}
  double GetY(){return Cart[1];}
  double GetTheta(){return Pol[0];}
  double GetR(){return Pol[1];}
  int GetName(){return Name;}

  void CarfromPol(){Cart[0]=Pol[1]*cos(Pol[0]); Cart[1]=Pol[1]*sin(Pol[0]);}
  
};

class Cromosome{
private:
  vector<Gene> Cstart;
  vector<Gene> C;
  double Prob;
  Random rnd;
public:
  Cromosome();
  Cromosome(Random,vector<Gene>);
  ~Cromosome();

  void SetProb(double x){Prob=x;}
  void SetGene(int i,Gene G){C[i]=G;}
  void NewC();
  double GetProb(){return Prob;}
  double Length();
  vector<Gene> GetC(){return C;}

  int Pbc(int);
};

class Population{
private:
  vector<Cromosome> P;
  vector<Cromosome> Pnew; 
  vector<double> LC,Prob;
  Random rnd;
  bool Ordinated=false;
public:
  Population();
  Population(Random r){rnd=r;};
  ~Population();

  void AddMember(Cromosome C){P.push_back(C);};
void ElimLastNElem(int i){for(int j=0;j<i;j++){P.pop_back();}};
  void CalcProb();
  void OrdinateAll();
	
  void Selection();
  void Crossover();
  void Mutation();
  void Evolve();

  double Lmean();
  int Pbc(int);
  vector<double> GetProb(){return Prob;}
  vector<double> GetLC(){return LC;}
  vector<Cromosome> GetP(){return P;}
  Cromosome GetBest();
  vector<Cromosome> GetNBest(int);
  Cromosome GetMember(int i){return P[i];}
  Cromosome GetMemberN(int i){return Pnew[i];}

};

class Error{
  private:

  protected:

  public:
    // constructors
    Error();
    // destructor
    ~Error();
    // methods
    double ErrStat(double, double, int);
    double Som(double a,double b);
    int add(int x, int y){
      return x+y;}
};



class FileInt{
  private:

  protected:

  public:
    // constructors
    FileInt();
    // destructor
    ~FileInt();
    // methods
    /*void WriteDouble(ofstream, double);*/
    void CreateCleanFile(string namefile );
};

class Operation{
  private:
    double mean=0,mean2=0, counter=0, counter2=0, min=0, max=0;
  protected:

  public:
    // constructors
    Operation();
    // destructor
    ~Operation();
    // methods
    void ResetMean();
    void UpdateMean(double x);
    double GetMean();
    void ResetMean2();
    void UpdateMean2(double x);
    double GetMean2();
    double ErrStat();
    void SetMinMax(double min,double max);
    void UpdateMinMax(double x);
    double GetMin();
    double GetMax();
};
#endif // __Library__


