#include "Edge.h"

#include <math.h>
#include<iostream>
#include<string>

Edge::Edge(){};
Edge::Edge(double b, unsigned int s1, unsigned int s2, Rcpp::String cstt, double param, double K, double a, double mini, double maxi) :
    beta(fabs(b)), state1(s1), state2(s2), constraint(cstt), parameter(fabs(param)), KK(K), aa(a), minn(mini), maxx(maxi){}


double Edge::getBeta() const {return(beta);}
unsigned int Edge::getState1() const {return(state1);}
unsigned int Edge::getState2() const {return(state2);}
std::string Edge::getConstraint() const {return(constraint);}
double Edge::getParameter() const {return(parameter);}
double Edge::getKK() const {return(KK);}
double Edge::getAA() const {return(aa);}
double Edge::getMinn() const {return(minn);}
double Edge::getMaxx() const {return(maxx);}

void Edge::show() const
{
  //std::cout << " state1: " << state1 << " state2: " << state2 << " constraint: " << constraint << " penalty: " << beta << " parameter: " << parameter << std::endl;
}
