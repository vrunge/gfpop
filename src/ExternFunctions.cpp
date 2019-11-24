#include "ExternFunctions.h"
#include <functional>
#include "math.h"

std::function<double*(Point const& pt)> cost_coeff;
std::function<double(const Cost&, double value)> cost_eval;

std::function<double(const Cost&)> cost_min;
std::function<double(const Cost&, Interval inter)> cost_minInterval;
std::function<double(const Cost&)> cost_argmin;
std::function<double(const Cost&, Interval inter)> cost_argminInterval;
std::function<double(const Cost&, Interval inter)> cost_argminBacktrack;

std::function<void(Cost& cost, double parameter)> cost_shift;
std::function<double(double bound, double parameter)> cost_interShift;
std::function<void(Cost& cost, double parameter)> cost_expDecay;
std::function<double(double bound, double parameter)> cost_interExpDecay;

std::function<Interval(const Cost&, double& level)> cost_intervalInterRoots;
std::function<int(const Cost&)> cost_age;
std::function<Interval()> cost_interval;
