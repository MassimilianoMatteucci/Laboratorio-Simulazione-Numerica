

#ifndef __Library__
#define __Library__
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
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

#endif // __Library__


