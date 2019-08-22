#ifndef DATA_H
#define DATA_H

#include "Rcpp.h"

#include <string>

struct Point
{
  double y;
  double w;
};

class Data
{
  public:
    Data();
    ~Data();

    void copy(Rcpp::NumericVector vectData, Rcpp::NumericVector vectWeight, unsigned int nd, unsigned int nw);

    double getm() const;
    double getM() const;
    unsigned int getn() const;

    double* gety() const;
    double* getw() const;
    Point* getVecPt() const;

    void show() const;

  private:
    Point* vecPt = NULL; ///see why we assign NULL in desctructor
    double m; ///min value
    double M; ///max value
    unsigned int n; ///data length
};

#endif // DATA_H
