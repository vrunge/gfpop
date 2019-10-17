#include "Data.h"

#include <string>
#include <iostream>
#include <fstream>

Data::Data(){}
Data::~Data()
{
  if (vecPt != NULL)
  {
    delete [] vecPt;
  }
}

//####### copy #######////####### copy #######////####### copy #######//
//####### copy #######////####### copy #######////####### copy #######//

void Data::copy(Rcpp::NumericVector vectData, Rcpp::NumericVector vectWeight, unsigned int nd, unsigned int nw)
{
  n = nd;
  vecPt = new Point[n]; ///array of Point of size n

  if(nw == n)
  {
    for(unsigned int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = vectWeight[i];}
  }
  else
  {
    for(unsigned int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = 1;}
  }

  /// Find the min and max values
  double mm = vectData[0];
  double MM = vectData[0];
  for(unsigned int i = 1 ; i < n ; i++)
  {
    if(vecPt[i].y < mm){mm = vecPt[i].y;}if(vecPt[i].y > MM){MM = vecPt[i].y;}
  }
  m = mm;
  M = MM;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//


double Data::getm() const {return(m);}
double Data::getM() const {return(M);}
unsigned int Data::getn() const {return(n);}

double* Data::gety() const
{
  double* y = new double[n];
  for (unsigned int i = 0 ; i < n ; i++){y[i] = vecPt[i].y;}
  return(y);
}

double* Data::getw() const
{
  double* w = new double[n];
  for (unsigned int i = 0 ; i < n ; i++){w[i] = vecPt[i].w;}
  return(w);
}

Point* Data::getVecPt() const {return(vecPt);}
