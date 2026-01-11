#ifndef EDGE_H
#define EDGE_H

#include<string>
#include "Rcpp.h"
#include "math.h" //to use INFINITY

class Edge
{
  public:
    Edge();
    Edge(unsigned int s1, unsigned int s2, Rcpp::String cstt = "std", double param = 0, double b = 0, double K = INFINITY, double a = 0, double mini = -INFINITY, double maxi = INFINITY);

    double getBeta() const;
    unsigned int getState1() const;
    unsigned int getState2() const;
    std::string getConstraint() const;
    double getParameter() const;
    double getKK() const;
    double getAA() const;
    double getMinn() const;
    double getMaxx() const;

    void show() const;
    bool isActive(unsigned int t) const;

  private:
    unsigned int state1;
    unsigned int state2;
    std::string constraint;
    double parameter;
    double beta;
    double KK;
    double aa;
    double minn;
    double maxx;
    unsigned int start_time;
    unsigned int end_time;
};

#endif // EDGE_H
