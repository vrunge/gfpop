#include "Data.h"
#include "termcolor.h"

#include <string>
#include <iostream>
#include <fstream>

Data::Data(){}
Data::Data(std::string file): file(file){}

Data::~Data()
{
  //std::cout << "Destructor called " << std::endl;
  if (vecPt != NULL)
  {
    delete [] vecPt;
  }
}

//####### copy #######////####### copy #######////####### copy #######//
//####### copy #######////####### copy #######////####### copy #######//


void Data::copy(Rcpp::NumericVector vectData, Rcpp::NumericVector vectWeight, int nd, int nw)
{
  n = nd;
  vecPt = new Point[n];

  if(nw == n)
  {
    for (int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = vectWeight[i];}
  }
  else
  {
    for (int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = 1;}
  }

  double mm = vecPt[0].y;
  double MM = vecPt[0].y;
  for (int i = 1 ; i < n ; i++)
  {
    if(vecPt[i].y < mm){mm = vecPt[i].y;}
    if(vecPt[i].y > MM){MM = vecPt[i].y;}
  }
  m = mm;
  M = MM;

}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//


double Data::getm() const {return(m);}
double Data::getM() const {return(M);}
int Data::getn() const {return(n);}



void  Data::setm(double mini){m = mini;}
void  Data::setM(double maxi){M = maxi;}


double* Data::gety() const
{
  double* y = new double[n];
  for (int i = 0 ; i < n ; i++){y[i] = vecPt[i].y;}
  return(y);
}

double* Data::getw() const
{
  double* w = new double[n];
  for (int i = 0 ; i < n ; i++){w[i] = vecPt[i].w;}
  return(w);
}

Point* Data::getVecPt() const {return(vecPt);}

//####### show #######////####### show #######////####### show #######//
//####### show #######////####### show #######////####### show #######//


void Data::show() const
{
  //std::cout << termcolor::red << "vector y" << termcolor::reset << std::endl;

  for (int i = 0 ; i < n ; i++)
  {
    //std::cout << vecPt[i].y << std::endl;
  }

  //std::cout << termcolor::red << "vector w" << termcolor::reset << std::endl;
  for (int i = 0 ; i < n ; i++)
  {
    //std::cout << vecPt[i].w << std::endl;
  }

  //std::cout << termcolor::red << "minimum : "<< m << termcolor::reset << std::endl;
  //std::cout << termcolor::red << "maximum : "<< M << termcolor::reset << std::endl;
  //std::cout << termcolor::red << "n : "<< n << termcolor::reset << std::endl;
}
