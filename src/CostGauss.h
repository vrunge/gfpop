#ifndef COSTGAUSS_H
#define COSTGAUSS_H

#include "Interval.h"
#include "Robust.h"

#include<vector>

struct Point
{
  double y;
  double w;
};

class CostGauss
{
  public:
    CostGauss();

    double getM_A() const;
    double getM_B() const;
    double getConstant() const;

    double minimum() const;
    double arg_minimum() const;
    double arg_maximum() const;

    double point_eval(double number) const;

    Interval intervalInterRoots() const;
    Interval intervalMinLess(double input, double bound, bool constPiece) const;

    int sign_Q2_Minus_Q1(CostGauss const& cost_Q2, double leftBound) const;

    ///OPERATIONS
    void operator+=(Point const& pt);
    void operator+=(double number);
    void shift(double parameter);
    void axisSymmetry();
    void opposition();

    CostGauss minus(CostGauss const& mycost);
    bool isEqual(CostGauss const& mycost) const;

    void addL1(Point const& pt, Robust const& robust, int slope);

    void show() const;
    void show(Interval const& inter) const;

    //Cost* operator-(Cost* cost1, Cost* cost2);

  private:
    double m_A;
    double m_B;
    double constant;
};


#endif // COSTGAUSS_H
