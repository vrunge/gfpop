#ifndef COST_H
#define COST_H

#include <functional>
#include"math.h"
#include"Interval.h"

struct Cost
{
  double m_A;
  double m_B;
  double constant;
  Cost();
  Cost(double* coeff);
};

void showCost(Cost& cost);
double log_factorial(double n);
double log_choose(double x, double n);

void addmyConstant(Cost& cost, double& cst);
void addCost(Cost& cost, const Cost& cost2);

std::function<double(const Cost&)> min_factory(const std::string& type);
std::function<double(const Cost&)> argmin_factory(const std::string& type);
std::function<double*(double)> coeff_factory(const std::string& type);
std::function<Interval(const Cost&, double& level)> intervalInterRoots_factory(const std::string& type);
std::function<int(const Cost&)> age_factory(const std::string& type);

#endif // COST_H
