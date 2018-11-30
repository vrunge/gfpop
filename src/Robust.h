#ifndef ROBUST_H
#define ROBUST_H

#include <string>
#include "math.h"

class Robust
{
  public:
    Robust(double K = INFINITY, double a = 0);

    double getThreshold() const;
    double getSlope() const;
    std::string getRobustType() const;

    void setThreshold(double K);
    void setSlope(double a);
    void findRobustType();

    void show() const;

  private:
    double threshold; ///Threshold (K) for outliers version of cost
    double slope; ///for theta greater that min + K. A L1 loss with slope > 0
    std::string robustType;

};

#endif // ROBUST_H
