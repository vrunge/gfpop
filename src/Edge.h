#ifndef EDGE_H
#define EDGE_H

#include<string>
#include "Rcpp.h"

class Edge
{
  public:
    Edge();
    Edge(double b, int s1 = 0, int s2 = 0, Rcpp::String cstt = "std", double param = 0);

    double getBeta() const;
    int getState1() const;
    int getState2() const;
    std::string getConstraint() const;
    double getParameter() const;

    void show() const;

  private:
    double beta;
    int state1;
    int state2;
    std::string constraint;
    double parameter;
};

#endif // EDGE_H
