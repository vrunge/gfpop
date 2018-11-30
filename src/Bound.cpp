#include "Bound.h"


Bound::Bound(double mini, double maxi, bool isConstr) : m(mini), M(maxi), isConstrained(isConstr){}

double Bound::getm() const {return(m);}
double Bound::getM() const {return(M);}
bool Bound::getIsConstrained() const {return(isConstrained);}

void Bound::setm(double mini){m = mini;}
void Bound::setM(double maxi){M = maxi;}
void Bound::setIsConstrained(bool isConstr){isConstrained = isConstr;}
