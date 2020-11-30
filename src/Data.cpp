#include "Data.h"

Data::Data(){}
Data::~Data()
{
  if (vecPt != NULL)
  {
    delete[] vecPt;
  }
}

//####### copy #######////####### copy #######////####### copy #######//
//####### copy #######////####### copy #######////####### copy #######//

void Data::copy(Rcpp::NumericVector vectData, Rcpp::NumericVector vectWeight, unsigned int nd, unsigned int nw)
{
  n = nd;
  vecPt = new Point[n]; ///array of Point of size n

  if(nw == n)
    {for(unsigned int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = vectWeight[i];}}
  else
    {for(unsigned int i = 0 ; i < n; i++){vecPt[i].y = vectData[i]; vecPt[i].w = 1;}}
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

unsigned int Data::getn() const {return(n);}
Point* Data::getVecPt() const {return(vecPt);}
