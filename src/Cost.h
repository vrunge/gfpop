#ifndef COST_H
#define COST_H

#include <functional>
#include "math.h"
#include "Interval.h"
#include "Data.h"


struct Cost
{
  double m_A;
  double m_B;
  double constant;
  Cost();
  Cost(double* coeff);
};

void addConstant(Cost& cost, double& cst);
void addCost(Cost& cost, const Cost& cost2);
Cost minusCost(Cost& cost, const Cost& cost2);
bool isEqual(Cost const& cost1, Cost const& cost2);

void showCost(const Cost& cost);

int signValue(double value);
double log_factorial(double n);
double log_choose(double x, double n);



///

std::function<double*(const Point&)> coeff_factory(const std::string& type);
std::function<double(const Cost&, double value)> eval_factory(const std::string& type);

std::function<double(const Cost&)> min_factory(const std::string& type);
std::function<double(const Cost&, Interval inter)> minInterval_factory(const std::string& type);
std::function<double(const Cost&)> argmin_factory(const std::string& type);

std::function<void(Cost& cost, double parameter)> shift_factory(const std::string& type);
std::function<double(double bound, double parameter)> interShift_factory(const std::string& type);
std::function<void(Cost& cost, double parameter)> expDecay_factory(const std::string& type);
std::function<double(double bound, double parameter)> interExpDecay_factory(const std::string& type);

std::function<Interval(const Cost&, double& level)> intervalInterRoots_factory(const std::string& type);

std::function<int(const Cost&)> age_factory(const std::string& type);
std::function<Interval()> interval_factory(const std::string& type);

#endif // COST_H
