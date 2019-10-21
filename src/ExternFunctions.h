#ifndef EXTERNFUNCTIONS_H
#define EXTERNFUNCTIONS_H

#include <functional>
#include"Cost.h"

extern std::function<double*(Point const& pt)> cost_coeff;
extern std::function<double(const Cost&)> cost_min;
extern std::function<double(const Cost&, Interval inter)> cost_minInterval;
extern std::function<double(const Cost&)> cost_argmin;
extern std::function<Interval(const Cost&, double& level)> cost_intervalInterRoots;
extern std::function<int(const Cost&)> cost_age;
extern std::function<Interval()> cost_interval;

#endif // EXTERNFUNCTIONS_H
