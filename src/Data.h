#ifndef DATA_H
#define DATA_H

#include "CostGauss.h"
#include "Rcpp.h"

#include <string>


class Data
{
  public:
    Data();
    Data(std::string file);
    ~Data();

    void read(bool weight = false);
    void copy(Rcpp::NumericVector vectData, Rcpp::NumericVector vectWeight, int nd, int nw);

    double getm() const;
    double getM() const;
    int getn() const;

    void setm(double mini);
    void setM(double maxi);


    double* gety() const;
    double* getw() const;
    Point* getVecPt() const;

    void show() const;

  private:
    std::string file;
    Point* vecPt = NULL; ///see why we assign NULL in desctructor
    double m;
    double M;
    int n;
};

#endif // DATA_H
