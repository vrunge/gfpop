#include "ExternFunctions.h"
#include <functional>
#include "Cost.h"

#include "math.h"

std::function<double*(Point const& pt)> cost_coeff;
std::function<double(const Cost&)> cost_min;
std::function<double(const Cost&, Interval inter)> cost_minInterval;
std::function<double(const Cost&)> cost_argmin;
std::function<Interval(const Cost&, double& level)> cost_intervalInterRoots;
std::function<int(const Cost&)> cost_age;
std::function<Interval()> cost_interval;
