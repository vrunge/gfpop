#include"Cost.h"

#include <vector>
#include <iostream>
#include <string>

Cost::Cost()
{
  m_A = 0;
  m_B = 0;
  constant = 0;
}

Cost::Cost(double* coeff)
{
  m_A = coeff[0];
  m_B = coeff[1];
  constant = coeff[2];
}

void addConstant(Cost& cost, double& cst){cost.constant = cost.constant + cst;}
Cost addCost(Cost& cost1, const Cost& cost2)
{
  Cost res = Cost();
  cost1.m_A = cost1.m_A + cost2.m_A;
  cost1.m_B = cost1.m_B + cost2.m_B;
  cost1.constant = cost1.constant + cost2.constant;
  return(res);
}

Cost minusCost(Cost& cost1, const Cost& cost2)
{
  Cost res = Cost();
  res.m_A = cost1.m_A - cost2.m_A;
  res.m_B = cost1.m_B - cost2.m_B;
  res.constant = cost1.constant - cost2.constant;
  return(res);
}

bool isEqual(Cost const& cost1, Cost const& cost2)
{
  return((cost1.m_A == cost2.m_A) && (cost1.m_B == cost2.m_B) && (cost1.constant == cost2.constant));
}

bool isConstant(Cost const& cost)
{
  return((cost.m_A == 0) && (cost.m_B == 0));
}

//####### showCost #######////####### showCost #######////####### showCost #######//
//####### showCost #######////####### showCost #######////####### showCost #######//

void showCost(const Cost& cost)
{
  //std::cout << " A: " << cost.m_A << " B: " << cost.m_B << " C: " << cost.constant << std::endl;
}

//####### simplefunctions #######////####### simplefunctions #######////####### simplefunctions #######//
//####### simplefunctions #######////####### simplefunctions #######////####### simplefunctions #######//

int signValue(double value)
{
  int res = 1;
  if(value < 0){res = -1;}
  return(res);
}

double log_factorial(double n)
{
  double res = 0;
  for(int i = 2; i < floor(n) + 1; i++){res = res + log((double)(i));}
  return(res);
}

double log_choose(double x, double n)
{
  if(n == x){return(0);}else if(x == 0 && n != 0){return(0);}
  else if(x == 1){return(log(n));}
  else
  {
    return(log_factorial(n) - log_factorial(x) - log_factorial(n - x));
  }
}

//////////////
//////////////
////////////// COEFF + EVAL
//////////////
//////////////

//####### coefficients #######////####### coefficients #######////####### coefficients #######//
//####### coefficients #######////####### coefficients #######////####### coefficients #######//

///mean cost = m_A*THETA^2 + m_B*THETA + constant
///variance cost = m_A*THETA - m_B*log(THETA) + constant
/// poisson cost = m_A*THETA - m_B*log(THETA) + constant
/// exp cost = m_A*THETA - m_B*log(THETA) + constant
/// negbin cost = - m_A*log(THETA) - m_B*log(1-THETA) + constant

double* mean_coeff(Point const& pt)
{
  double* coeff = new double[3];
  coeff[0] = pt.w;
  coeff[1] = -2 * pt.w * pt.y;
  coeff[2] = pt.w * pt.y * pt.y;
  return(coeff);
}

double* variance_coeff(Point const& pt)
{
  double* coeff = new double[3];
  coeff[0] = pt.w * pt.y * pt.y;
  coeff[1] = pt.w;
  //coeff[2] = log(2 * M_PI);
  coeff[2] = 0;
  return(coeff);
}

double* poisson_coeff(Point const& pt)
{
  double* coeff = new double[3];
  coeff[0] = pt.w;
  coeff[1] = pt.w * pt.y;
  //coeff[2] = log_factorial(y);
  coeff[2] = 0;
  return(coeff);
}

double* exp_coeff(Point const& pt)
{
  double* coeff = new double[3];
  coeff[0] = pt.w * pt.y;
  coeff[1] = pt.w;
  coeff[2] = 0;
  return(coeff);
}

double* negbin_coeff(Point const& pt)
{
  double* coeff = new double[3];
  coeff[0] = pt.w;
  coeff[1] = pt.w * pt.y;
  coeff[2] = 0;
  return(coeff);
}

//####### eval #######////####### eval #######////####### eval #######//
//####### eval #######////####### eval #######////####### eval #######//
/// NEVER = - INFINITY by construction

double mean_eval(const Cost& cost, double value)
{
  double res = INFINITY; ///default case
  if(value != -INFINITY && value != INFINITY){res = cost.m_A * value * value + cost.m_B * value + cost.constant;}
  else if(cost.m_A == 0 && cost.m_B == 0){res = cost.constant;}
  return(res);
}

double variance_eval(const Cost& cost, double value)
{
  double res = INFINITY; ///default case
  if(value != 0 && value != INFINITY){res = cost.m_A * value - cost.m_B * log(value) + cost.constant;}
  else if(value == 0 && cost.m_B == 0){res = cost.constant;}
  else if(cost.m_A == 0 && cost.m_B == 0){res = cost.constant;}
  return(res);
}

double negbin_eval(const Cost& cost, double value)
{
  double res = INFINITY; ///default case
  if(value != 0 && value != 1){res = - cost.m_A * log(value) - cost.m_B * log(1 - value) + cost.constant;}
  else if(value == 0 && cost.m_A == 0){res = cost.constant;}
  else if(value  == 1 && cost.m_B == 0){res = cost.constant;}
  else if(cost.m_A == 0 && cost.m_B == 0){res = cost.constant;}
  return(res);
}

//////////////
//////////////
////////////// MIN + ARGMIN
//////////////
//////////////

//####### minimum #######////####### minimum #######////####### minimum #######//
//####### minimum #######////####### minimum #######////####### minimum #######//

double mean_min(const Cost& cost)
{
  double res = -INFINITY;
  if(cost.m_A > 0){res = - (cost.m_B * cost.m_B/(4 * cost.m_A)) + cost.constant;}
  else if(cost.m_A == 0 && cost.m_B == 0){res = cost.constant;}
  return(res);
}

double variance_min(const Cost& cost)
{
  double res = -INFINITY;
  if(cost.m_A > 0 && cost.m_B > 0){res = cost.m_B - cost.m_B * log(cost.m_B/cost.m_A) + cost.constant;}
  else if(cost.m_A >= 0 && cost.m_B == 0){res = cost.constant;}
  return(res);
}

double negbin_min(const Cost& cost)
{
  double res = -INFINITY;
  if(cost.m_A > 0 && cost.m_B > 0){res = - cost.m_A * log(cost.m_A/(cost.m_A + cost.m_B)) - cost.m_B * log(cost.m_B/(cost.m_A + cost.m_B)) + cost.constant;}
  else if(cost.m_A == 0 || cost.m_B == 0){res = cost.constant;}
  return(res);
}


//####### minInterval #######////####### minInterval #######////####### minInterval #######//
//####### minInterval #######////####### minInterval #######////####### minInterval #######//

double mean_minInterval(const Cost& cost, Interval inter)
{
  double minimum = -INFINITY;
  //case m_A > 0
  if(cost.m_A > 0)
  {
    minimum = - (cost.m_B * cost.m_B/(4 * cost.m_A)) + cost.constant;
    double argmin = - cost.m_B/(2 * cost.m_A);
    if(argmin < inter.geta()){minimum = cost.m_A * inter.geta() * inter.geta() + cost.m_B * inter.geta() + cost.constant;}
    if(argmin > inter.getb()){minimum = cost.m_A * inter.getb() * inter.getb() + cost.m_B * inter.getb() + cost.constant;}
  }

  //case m_A = 0 & m_B != 0
  else if((cost.m_A == 0) && (cost.m_B != 0))
  {
    if(cost.m_B > 0){minimum = cost.m_B * inter.geta() + cost.constant;}
      else{minimum = cost.m_B * inter.getb() + cost.constant;}
  }
  //case m_A = 0 & m_B = 0
  else if((cost.m_A == 0) && (cost.m_B == 0)){minimum = cost.constant;}

  return(minimum);
}


double variance_minInterval(const Cost& cost, Interval inter)
{
  double minimum = -INFINITY;
  //case m_A > 0 & cost.m_B > 0
  if(cost.m_A > 0 && cost.m_B > 0)
  {
    minimum = cost.m_B - cost.m_B * log(cost.m_B/cost.m_A) + cost.constant;
    double argmin = cost.m_B/cost.m_A;
    if(argmin < inter.geta()){minimum = cost.m_A * inter.geta() - cost.m_B * log(inter.geta()) + cost.constant;}
    if(argmin > inter.getb()){minimum = cost.m_A * inter.getb() - cost.m_B * log(inter.getb()) + cost.constant;}
  }
  //case m_A != 0 & m_B = 0
  else if((cost.m_A != 0) && (cost.m_B == 0))
  {
    if(cost.m_A > 0){minimum = cost.m_A * inter.geta() + cost.constant;}
      else{minimum = cost.m_A * inter.getb() + cost.constant;}
  }
  //case m_A = 0 & m_B = 0
  else if((cost.m_A == 0) && (cost.m_B == 0)){minimum = cost.constant;}

  return(minimum);
}


double negbin_minInterval(const Cost& cost, Interval inter)
{
  double minimum = -INFINITY;
  //case m_A > 0 & cost.m_B > 0
  if(cost.m_A > 0 && cost.m_B > 0)
  {
    minimum = - cost.m_A * log(cost.m_A/(cost.m_A + cost.m_B)) - cost.m_B * log(cost.m_B/(cost.m_A + cost.m_B)) + cost.constant;
    double argmin = cost.m_A/(cost.m_A + cost.m_B);
    if(argmin < inter.geta()){minimum = - cost.m_A * log(inter.geta()) - cost.m_B * log(1 - inter.geta()) + cost.constant;}
    if(argmin > inter.getb()){minimum = - cost.m_A * log(inter.getb()) - cost.m_B * log(1 - inter.getb()) + cost.constant;}
  }
  else if(cost.m_A > 0 && cost.m_B == 0){minimum = - cost.m_A * log(inter.getb()) + cost.constant;}
  else if(cost.m_A == 0 && cost.m_B > 0){minimum = - cost.m_B * log(1 - inter.geta()) + cost.constant;}
  else if((cost.m_A == 0) && (cost.m_B == 0)){minimum = cost.constant;}

  return(minimum);
}

//####### argmin #######////####### argmin #######////####### argmin #######//
//####### argmin #######////####### argmin #######////####### argmin #######//

double mean_argmin(const Cost& cost)
{
  double argmin = INFINITY;
  if(cost.m_A == 0)
  {
    //if(m_B<=0){argmin = INFINITY;}
    if(cost.m_B > 0){argmin = -INFINITY;}
  }
  else
  {
    argmin = - cost.m_B/(2 * cost.m_A);
  }
  return(argmin);
}

double variance_argmin(const Cost& cost)
{
  double argmin = INFINITY;
  if(cost.m_B == 0)
  {
    //if(m_A <= 0){argmin = INFINITY;}
    if(cost.m_A > 0){argmin = 0;}
  }
  else
  {
    argmin = cost.m_B/cost.m_A;
  }
  return(1/argmin); ///for the variance parameter
}


double poisson_argmin(const Cost& cost)
{
  double argmin = INFINITY;
  if(cost.m_B == 0)
  {
    //if(m_A <= 0){argmin = INFINITY;}
    if(cost.m_A > 0){argmin = 0;}
  }
  else if(cost.m_A != 0)
  {
    argmin = cost.m_B/cost.m_A; ///for the functional cost parameter
  }
  return(argmin);
}


double negbin_argmin(const Cost& cost)
{
  double argmin = 0.5;
  if(cost.m_A > 0 && cost.m_B > 0){argmin = cost.m_A/(cost.m_A + cost.m_B);}
  else if(cost.m_A > 0 && cost.m_B == 0){argmin = 1;}
  else if(cost.m_A == 0 && cost.m_B > 0){argmin = 0;}
  return(argmin);
}

//####### argminInterval #######////####### argminInterval #######////####### argminInterval #######//
//####### argminInterval #######////####### argminInterval #######////####### argminInterval #######//

double mean_argminInterval(const Cost& cost, Interval inter)
{
  double argmin = inter.getb();
  if(cost.m_A == 0)
  {
    if(cost.m_B == 0){argmin = (inter.geta() + inter.getb())/2;}
    else if(cost.m_B > 0){argmin = inter.geta();}
  }
  else
  {
    argmin = - cost.m_B/(2 * cost.m_A);
    if(argmin < inter.geta()){argmin = inter.geta();}
    else if(argmin > inter.getb()){argmin = inter.getb();}
  }
  return(argmin);
}

double variance_argminInterval(const Cost& cost, Interval inter)
{
  double argmin = inter.getb();
  if(cost.m_B == 0)
  {
    if(cost.m_A == 0){argmin = (inter.geta() + inter.getb())/2;}
    else if(cost.m_A > 0){argmin = inter.geta();}
  }
  else
  {
    argmin = cost.m_B/cost.m_A;
    if(argmin < inter.geta()){argmin = inter.geta();}
    else if(argmin > inter.getb()){argmin = inter.getb();}
  }
  return(1/argmin); ///for the variance parameter
}


double poisson_argminInterval(const Cost& cost, Interval inter)
{
  double argmin = inter.getb();
  if(cost.m_B == 0)
  {
    if(cost.m_A == 0){argmin = (inter.geta() + inter.getb())/2;}
    else if(cost.m_A > 0){argmin = inter.geta();}
  }
  else
  {
    argmin = cost.m_B/cost.m_A;
    if(argmin < inter.geta()){argmin = inter.geta();}
    else if(argmin > inter.getb()){argmin = inter.getb();}
  }
  return(argmin);
}


double negbin_argminInterval(const Cost& cost, Interval inter)
{
  double argmin = inter.getb();
  if(cost.m_A > 0 && cost.m_B > 0)
  {
    argmin = cost.m_A/(cost.m_A + cost.m_B);
    if(argmin < inter.geta()){argmin = inter.geta();}
    else if(argmin > inter.getb()){argmin = inter.getb();}
  }

  else if(cost.m_A > 0 && cost.m_B == 0){argmin = inter.getb();}
  else if(cost.m_A == 0 && cost.m_B > 0){argmin = inter.geta();}
  else if(cost.m_A == 0 && cost.m_B == 0){argmin = (inter.geta() + inter.getb())/2;}
  return(argmin);
}
//////////////
//////////////
////////////// TRANSFORMATIONS
//////////////
//////////////

//####### shift #######////####### shift #######////####### shift #######//
//####### shift #######////####### shift #######////####### shift #######//

void mean_shift(Cost& cost, double parameter)
{
  cost.constant = cost.constant + parameter * (cost.m_A * parameter - cost.m_B);
  cost.m_B = cost.m_B - 2 * cost.m_A * parameter;
}

void variance_shift(Cost& cost, double parameter)
{
  if(parameter > 0)
  {
    cost.m_A = cost.m_A / parameter;
    cost.constant = cost.constant + cost.m_B * log(parameter);
  }
  else if(parameter < 0)
  {
    cost.m_A = cost.m_A * fabs(parameter);
    cost.constant = cost.constant - cost.m_B * log(fabs(parameter));
  }
}

void negbin_shift(Cost& cost, double parameter){}


//####### interShift #######////####### interShift #######////####### interShift #######//
//####### interShift #######////####### interShift #######////####### interShift #######//

double mean_interShift(double bound, double parameter)
{
  return(bound + parameter);
}

double variance_interShift(double bound, double parameter)
{
  double res = bound;
  if(parameter > 0){res = bound * parameter;}
  else if(parameter < 0){res = bound / fabs(parameter);}
  return(res);
}

double negbin_interShift(double bound, double parameter){return(bound);}



//####### expDecay #######////####### expDecay #######////####### expDecay #######//
//####### expDecay #######////####### expDecay #######////####### expDecay #######//

void mean_expDecay(Cost& cost, double gamma)
{
  cost.m_A = cost.m_A / (gamma * gamma);
  cost.m_B = cost.m_B / gamma;
}

void variance_expDecay(Cost& cost, double gamma)
{
  cost.m_A = cost.m_A / gamma;
  cost.constant = cost.constant + cost.m_B * log(gamma);
}

void negbin_expDecay(Cost& cost, double gamma){}



//####### interExpDecay #######////####### interExpDecay #######////####### interExpDecay #######//
//####### interExpDecay #######////####### interExpDecay #######////####### interExpDecay #######//

double mean_interExpDecay(double bound, double gamma)
{
  return(bound * gamma);
}

double variance_interExpDecay(double bound, double gamma)
{
  return(bound * gamma);
}

double negbin_interExpDecay(double bound, double gamma){return(bound);}


//####### intervalInterRoots #######////####### intervalInterRoots #######////####### intervalInterRoots #######//
//####### intervalInterRoots #######////####### intervalInterRoots #######////####### intervalInterRoots #######//

//####### MEAN #######////
Interval mean_intervalInterRoots(const Cost& cost, double& level)
{
  Interval newElement = Interval();
  double Delta = cost.m_B * cost.m_B - 4 * cost.m_A * (cost.constant - level);

  if(Delta > 0)
  {
    double R = sqrt(Delta);
    if(cost.m_A > 0){newElement = Interval((- cost.m_B - R)/(2 * cost.m_A), (- cost.m_B + R)/(2 * cost.m_A));}
    if(cost.m_A < 0){newElement = Interval((- cost.m_B + R)/(2 * cost.m_A), (- cost.m_B - R)/(2 * cost.m_A));}
  }
  if(cost.m_A == 0)
  {
    if(cost.m_B > 0){newElement = Interval(-INFINITY, -cost.constant/cost.m_B);}
    if(cost.m_B < 0){newElement = Interval(-cost.constant/cost.m_B, INFINITY);}
  }

  return(newElement);
}

//####### VARIANCE #######////
Interval variance_intervalInterRoots(const Cost& cost, double& level)
{
  Interval newElement = Interval();

  double eps = 0.000001;
  double temp = 1;
  //roots of A THETA - B log THETA + C = 0
  // <=> roots of THETA - log THETA = 1 + a
  double U = cost.m_A/cost.m_B;
  double a = -(((cost.constant - level)/cost.m_B) + log(U) + 1);

  if(a > 0)
  {
    int nb = 0;

    double leftRoot = -(1 + a);
    while(fabs(temp - leftRoot) > eps && nb < 100)
    {
      temp = leftRoot;
      leftRoot = leftRoot - 1 - (leftRoot + a)/(1 - exp(leftRoot));
      nb = nb + 1;
    }

    nb = 0;
    temp = 1;
    double rightRoot = 1 + a;
    while(fabs(temp - rightRoot) > eps && nb < 100)
    {
      temp = rightRoot;
      rightRoot = (log(rightRoot) + a)*rightRoot/(rightRoot - 1);
      nb = nb + 1;
    }

    leftRoot = exp(leftRoot)/U;
    rightRoot = rightRoot/U;

    newElement.seta(leftRoot);
    newElement.setb(rightRoot);
    if(leftRoot >= rightRoot){newElement = Interval();}
    }

  return(newElement);
}

//####### POISSON #######////
Interval poisson_intervalInterRoots(const Cost& cost, double& level)
{
  Interval newElement = Interval();

  if(cost.m_B > 0)
  {
    double eps = 0.000001;
    double temp = 1;
    //roots of A THETA - B log THETA + C = 0
    // <=> roots of THETA - log THETA = 1 + a
    double U = cost.m_A/cost.m_B;
    double a = -(((cost.constant - level)/cost.m_B) + log(U) + 1);

    if(a > 0)
    {
      int nb = 0;

      double leftRoot = -(1 + a);
      while(fabs(temp - leftRoot) > eps && nb < 100)
      {
        temp = leftRoot;
        leftRoot = leftRoot - 1 - (leftRoot + a)/(1 - exp(leftRoot));
        nb = nb + 1;
      }

      nb = 0;
      temp = 1;

      double rightRoot = 1 + a;
      while(fabs(temp - rightRoot) > eps && nb < 100)
      {
        temp = rightRoot;
        rightRoot = (log(rightRoot) + a)*rightRoot/(rightRoot - 1);
        nb = nb + 1;
      }

      leftRoot = exp(leftRoot)/U;
      rightRoot = rightRoot/U;

      newElement.seta(leftRoot);
      newElement.setb(rightRoot);
      if(leftRoot >= rightRoot){newElement = Interval();}
    }
  }
  else
  {
    newElement.seta(0);
    newElement.setb(level - cost.constant);
  }

  return(newElement);
}


//####### negbin #######////
Interval negbin_intervalInterRoots(const Cost& cost, double& level)
{
  Interval newElement = Interval();

  double eps = 0.000001;
  double temp = 1;
  //roots of - A log(THETA) - B log(1-THETA) + C = level
  double U = cost.m_A/(cost.m_A + cost.m_B);
  double a = level + cost.m_A*log(U) + cost.m_B*(1-U) - cost.constant;

  if(a > 0)
  {
    int nb = 0;

    double leftRoot = (cost.constant - level)/cost.m_A;
    while(fabs(temp - leftRoot) > eps && nb < 100)
    {
      temp = leftRoot;
      leftRoot = leftRoot - ((1 + exp(leftRoot))/(- cost.m_A + cost.m_B * exp(leftRoot))) * (- cost.m_A * leftRoot + (cost.m_A + cost.m_B)*log(1 + exp(leftRoot)) + cost.constant - level);
      nb = nb + 1;
    }

    nb = 0;
    temp = 1;

    double rightRoot = (level - cost.constant)/cost.m_B;
    while(fabs(temp - rightRoot) > eps && nb < 100)
    {
      temp = rightRoot;
      rightRoot = rightRoot - ((1 + exp(rightRoot))/(- cost.m_A + cost.m_B * exp(rightRoot))) * (- cost.m_A * rightRoot + (cost.m_A + cost.m_B)*log(1 + exp(rightRoot)) + cost.constant - level);
      nb = nb + 1;
    }

    leftRoot = exp(leftRoot)/(1 + exp(leftRoot));
    rightRoot = exp(rightRoot)/(1 + exp(rightRoot));

    newElement.seta(leftRoot);
    newElement.setb(rightRoot);
    if(leftRoot >= rightRoot){newElement = Interval();}
  }

  return(newElement);
}


//####### ages #######////####### ages #######////####### ages #######//
//####### ages #######////####### ages #######////####### ages #######//

int mean_age(const Cost& cost){return((int) cost.m_A);}
int variance_age(const Cost& cost){return((int) cost.m_B);}
int poisson_age(const Cost& cost){return((int) cost.m_A);}
int exp_age(const Cost& cost){return((int) cost.m_B);}
int negbin_age(const Cost& cost){return((int) cost.m_A);}


//####### intervals #######////####### intervals #######////####### intervals #######//
//####### intervals #######////####### intervals #######////####### intervals #######//

Interval mean_interval(){return(Interval(-INFINITY,INFINITY));}
Interval variance_interval(){return(Interval(0,INFINITY));}
Interval negbin_interval(){return(Interval(0,1));}


//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#
//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#
//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#
//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#
//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#//#

//####### factories #######////####### factories #######////####### factories #######//
//####### factories #######////####### factories #######////####### factories #######//

std::function<double*(const Point&)> coeff_factory(const std::string& type)
{
  std::function<double*(const Point&)> fct;
  if(type == "mean"){fct = std::function<double*(const Point&)>(mean_coeff);}
  if(type == "variance"){fct = std::function<double*(const Point&)>(variance_coeff);}
  if(type == "poisson"){fct = std::function<double*(const Point&)>(poisson_coeff);}
  if(type == "exp"){fct = std::function<double*(const Point&)>(exp_coeff);}
  if(type == "negbin"){fct = std::function<double*(const Point&)>(negbin_coeff);}
  return(fct);
}

std::function<double(const Cost&, double value)> eval_factory(const std::string& type)
{
  std::function<double(const Cost&, double value)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&, double value)>(mean_eval);}
  if(type == "variance"){fct = std::function<double(const Cost&, double value)>(variance_eval);}
  if(type == "poisson"){fct = std::function<double(const Cost&, double value)>(variance_eval);}
  if(type == "exp"){fct = std::function<double(const Cost&, double value)>(variance_eval);}
  if(type == "negbin"){fct = std::function<double(const Cost&, double value)>(negbin_eval);}
  return(fct);
}


//////////////
//////////////
////////////// MIN + ARGMIN
//////////////
//////////////

std::function<double(const Cost&)> min_factory(const std::string& type)
{
  std::function<double(const Cost&)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&)>(mean_min);}
  if(type == "variance"){fct = std::function<double(const Cost&)>(variance_min);}
  if(type == "poisson"){fct = std::function<double(const Cost&)>(variance_min);}
  if(type == "exp"){fct = std::function<double(const Cost&)>(variance_min);}
  if(type == "negbin"){fct = std::function<double(const Cost&)>(negbin_min);}
  return(fct);
}

std::function<double(const Cost&, Interval inter)> minInterval_factory(const std::string& type)
{
  std::function<double(const Cost&, Interval inter)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&, Interval inter)>(mean_minInterval);}
  if(type == "variance"){fct = std::function<double(const Cost&, Interval inter)>(variance_minInterval);}
  if(type == "poisson"){fct = std::function<double(const Cost&, Interval inter)>(variance_minInterval);}
  if(type == "exp"){fct = std::function<double(const Cost&, Interval inter)>(variance_minInterval);}
  if(type == "negbin"){fct = std::function<double(const Cost&, Interval inter)>(negbin_minInterval);}
  return(fct);
}


std::function<double(const Cost&)> argmin_factory(const std::string& type)
{
  std::function<double(const Cost&)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&)>(mean_argmin);}
  if(type == "variance"){fct = std::function<double(const Cost&)>(poisson_argmin);}
  if(type == "poisson"){fct = std::function<double(const Cost&)>(poisson_argmin);}
  if(type == "exp"){fct = std::function<double(const Cost&)>(poisson_argmin);}
  if(type == "negbin"){fct = std::function<double(const Cost&)>(negbin_argmin);}
  return(fct);
}


std::function<double(const Cost&, Interval inter)> argminInterval_factory(const std::string& type)
{
  std::function<double(const Cost&, Interval inter)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&, Interval inter)>(mean_argminInterval);}
  if(type == "variance"){fct = std::function<double(const Cost&, Interval inter)>(poisson_argminInterval);}
  if(type == "poisson"){fct = std::function<double(const Cost&, Interval inter)>(poisson_argminInterval);}
  if(type == "exp"){fct = std::function<double(const Cost&, Interval inter)>(poisson_argminInterval);}
  if(type == "negbin"){fct = std::function<double(const Cost&, Interval inter)>(negbin_argminInterval);}
  return(fct);
}

std::function<double(const Cost&, Interval inter)> argminBacktrack_factory(const std::string& type)
{
  std::function<double(const Cost&, Interval inter)> fct;
  if(type == "mean"){fct = std::function<double(const Cost&, Interval inter)>(mean_argminInterval);}
  if(type == "variance"){fct = std::function<double(const Cost&, Interval inter)>(variance_argminInterval);}
  if(type == "poisson"){fct = std::function<double(const Cost&, Interval inter)>(poisson_argminInterval);}
  if(type == "exp"){fct = std::function<double(const Cost&, Interval inter)>(poisson_argminInterval);}
  if(type == "negbin"){fct = std::function<double(const Cost&, Interval inter)>(negbin_argminInterval);}
  return(fct);
}

//////////////
//////////////
////////////// TRANSFORMATIONS
//////////////
//////////////

std::function<void(Cost& cost, double parameter)> shift_factory(const std::string& type)
{
  std::function<void(Cost& cost, double parameter)> fct;
  if(type == "mean"){fct = std::function<void(Cost& cost, double parameter)>(mean_shift);}
  if(type == "variance"){fct = std::function<void(Cost& cost, double parameter)>(variance_shift);}
  if(type == "poisson"){fct = std::function<void(Cost& cost, double parameter)>(variance_shift);}
  if(type == "exp"){fct = std::function<void(Cost& cost, double parameter)>(variance_shift);}
  if(type == "negbin"){fct = std::function<void(Cost& cost, double parameter)>(negbin_shift);}
  return(fct);
}


std::function<double(double bound, double parameter)> interShift_factory(const std::string& type)
{
  std::function<double(double bound, double parameter)> fct;
  if(type == "mean"){fct = std::function<double(double bound, double parameter)>(mean_interShift);}
  if(type == "variance"){fct = std::function<double(double bound, double parameter)>(variance_interShift);}
  if(type == "poisson"){fct = std::function<double(double bound, double parameter)>(variance_interShift);}
  if(type == "exp"){fct = std::function<double(double bound, double parameter)>(variance_interShift);}
  if(type == "negbin"){fct = std::function<double(double bound, double parameter)>(negbin_interShift);}
  return(fct);
}



std::function<void(Cost& cost, double gamma)> expDecay_factory(const std::string& type)
{
  std::function<void(Cost& cost, double gamma)> fct;
  if(type == "mean"){fct = std::function<void(Cost& cost, double gamma)>(mean_expDecay);}
  if(type == "variance"){fct = std::function<void(Cost& cost, double gamma)>(variance_expDecay);}
  if(type == "poisson"){fct = std::function<void(Cost& cost, double gamma)>(variance_expDecay);}
  if(type == "exp"){fct = std::function<void(Cost& cost, double gamma)>(variance_expDecay);}
  if(type == "negbin"){fct = std::function<void(Cost& cost, double gamma)>(negbin_expDecay);}
  return(fct);
}


std::function<double(double bound, double gamma)> interExpDecay_factory(const std::string& type)
{
  std::function<double(double bound, double gamma)> fct;
  if(type == "mean"){fct = std::function<double(double bound, double gamma)>(mean_interExpDecay);}
  if(type == "variance"){fct = std::function<double(double bound, double gamma)>(variance_interExpDecay);}
  if(type == "poisson"){fct = std::function<double(double bound, double gamma)>(variance_interExpDecay);}
  if(type == "exp"){fct = std::function<double(double bound, double gamma)>(variance_interExpDecay);}
  if(type == "negbin"){fct = std::function<double(double bound, double gamma)>(negbin_interExpDecay);}
  return(fct);
}


std::function<Interval(const Cost&, double& level)> intervalInterRoots_factory(const std::string& type)
{
  std::function<Interval(const Cost&, double& level)> fct;
  if(type == "mean"){fct = std::function<Interval(const Cost&, double& level)>(mean_intervalInterRoots);}
  if(type == "variance"){fct = std::function<Interval(const Cost&, double& level)>(variance_intervalInterRoots);}
  if(type == "poisson"){fct = std::function<Interval(const Cost&, double& level)>(poisson_intervalInterRoots);}
  if(type == "exp"){fct = std::function<Interval(const Cost&, double& level)>(variance_intervalInterRoots);}
  if(type == "negbin"){fct = std::function<Interval(const Cost&, double& level)>(negbin_intervalInterRoots);}
  return(fct);
}

std::function<int(const Cost&)> age_factory(const std::string& type)
{
  std::function<int(const Cost&)> fct;
  if(type == "mean"){fct = std::function<int(const Cost&)>(mean_age);}
  if(type == "variance"){fct = std::function<int(const Cost&)>(variance_age);}
  if(type == "poisson"){fct = std::function<int(const Cost&)>(poisson_age);}
  if(type == "exp"){fct = std::function<int(const Cost&)>(exp_age);}
  if(type == "negbin"){fct = std::function<int(const Cost&)>(negbin_age);}
  return(fct);
}


std::function<Interval()> interval_factory(const std::string& type)
{
  std::function<Interval()> fct;
  if(type == "mean"){fct = std::function<Interval()>(mean_interval);}
  if(type == "variance"){fct = std::function<Interval()>(variance_interval);}
  if(type == "poisson"){fct = std::function<Interval()>(variance_interval);}
  if(type == "exp"){fct = std::function<Interval()>(variance_interval);}
  if(type == "negbin"){fct = std::function<Interval()>(negbin_interval);}
  return(fct);
}

