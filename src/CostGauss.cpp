#include "CostGauss.h"

#include"math.h"
#include <iostream>
#include<vector>


CostGauss::CostGauss(): m_A(0), m_B(0), constant(0){}

double CostGauss::getM_A() const {return(m_A);}
double CostGauss::getM_B() const {return(m_B);}
double CostGauss::getConstant() const {return(constant);}


//####### min arg_min #######////####### min arg_min #######////####### min arg_min #######//
//####### min arg_min #######////####### min arg_min #######////####### min arg_min #######//

double CostGauss::minimum() const
{
  double minimum = -INFINITY;
  if(m_A > 0){minimum = -(m_B*m_B/(4*m_A)) + constant;}
  if((m_A == 0) && (m_B == 0)){minimum = constant;}

  return(minimum);
}

double CostGauss::arg_minimum() const
{
  double argmin = INFINITY;
  if(m_A == 0)
  {
    ///if(m_B<=0){argmin = INFINITY;}
    if(m_B > 0){argmin = -INFINITY;}
  }
  else
  {
    argmin = -m_B/(2*m_A);
  }

  return(argmin);
}
/// ARGMIN of a constant = + INFINITY


double CostGauss::arg_maximum() const
{
  double argmax = INFINITY;
  if(m_A == 0)
  {
    if(m_B < 0){argmax = -INFINITY;}
    ///if(m_B > 0){argmin = -INFINITY;}
  }
  else
  {
    argmax = -m_B/(2*m_A);
  }

  return(argmax);
}

//####### point_eval #######////####### point_eval #######////####### point_eval #######//
//####### point_eval #######////####### point_eval #######////####### point_eval #######//

double CostGauss::point_eval(double number) const {return((m_A*number*number) + (m_B*number) + constant);}


//####### comparisons #######////####### comparisons #######////####### comparisons #######//
//####### comparisons #######////####### comparisons #######////####### comparisons #######//

Interval CostGauss::intervalInterRoots() const
{
  Interval newElement = Interval();
  double Delta = m_B*m_B - 4*m_A*constant;

  if((m_A != 0) && (Delta >= 0))
  {
    double argmin = arg_minimum();
    double R = sqrt(Delta)/(2*m_A);
    newElement = Interval(argmin - fabs(R), argmin + fabs(R));
  }

  if(m_A == 0)
  {
    if((m_B == 0)&&(constant == 0)){newElement = Interval(-INFINITY, INFINITY);} ///!!!!cut = constant gives the whole real line
    if(m_B > 0){newElement = Interval(-INFINITY, -constant/m_B);}
    if(m_B < 0){newElement = Interval(-constant/m_B, INFINITY);}
  }

  return(newElement);
}


Interval CostGauss::intervalMinLess(double current_min, double bound, bool constPiece) const
{
  Interval response = Interval(); /// Interval = (INFINITY, INFINITY)
  double mini = minimum();

  if(current_min > mini)
  {
    double argmini = arg_minimum();
    if(bound < argmini)
    {
      if(constPiece == true)
      {
        CostGauss newCost = *this;
        newCost += -current_min;
        response = newCost.intervalInterRoots();
        response.setb(argmini);
      }
      else
      { /// i.e. point_eval(bound) == current_min : continuity condition
        response.seta(bound);
        response.setb(argmini);
      }
    }
  }

  return(response);
}


//####### sign_Q2_Minus_Q1 #######////####### sign_Q2_Minus_Q1 #######////####### sign_Q2_Minus_Q1 #######//
//####### sign_Q2_Minus_Q1 #######////####### sign_Q2_Minus_Q1 #######////####### sign_Q2_Minus_Q1 #######//


int CostGauss::sign_Q2_Minus_Q1(CostGauss const& cost_Q2, double leftBound) const
{
  int Q2_Minus_Q1 = 0;
  ///0 = positive sign

  double minQ1 = this -> point_eval(leftBound);
  double minQ2 = cost_Q2.point_eval(leftBound);

  if(minQ2 > minQ1){Q2_Minus_Q1 = 1;}
  if(minQ2 < minQ1){Q2_Minus_Q1 = -1;}

  if(minQ1 == minQ2)
  {
    double m_A2 = cost_Q2.getM_A();
    double m_B2 = cost_Q2.getM_B();

    if(2 * m_A * leftBound + m_B == 2 * m_A2 * leftBound + m_B2)
    {
      if(m_A <= m_A2){Q2_Minus_Q1 = 1;}else{Q2_Minus_Q1 = -1;}   /// m_A = m_A2 <=> priority to Q1 (the current cost, not the competitive new changepoint cost)
    }
    else
    {
      if(2 * m_A * leftBound + m_B < 2 * m_A2 * leftBound + m_B2){Q2_Minus_Q1 = 1;}
      else{Q2_Minus_Q1 = -1;}
    }
  }

  //Q2_Minus_Q1 = 1 OR -1
  return(Q2_Minus_Q1);

}





//####### operations #######////####### operations #######////####### operations #######//
//####### operations #######////####### operations #######////####### operations #######//

/// m_A*THETA^2 + m_B*THETA + constant
/// pt.w*(THETA^2 - 2*pt.y*THETA + pt.y*pt.y)

void CostGauss::operator+=(Point const& pt)
{
  m_A = m_A + pt.w;
  m_B = m_B - 2 * pt.y * pt.w;
  constant = constant + pt.w * pt.y * pt.y;
}

void CostGauss::operator+=(double number)
{constant = constant + number;}


void CostGauss::shift(double parameter)
{
  constant = constant + parameter * (m_A * parameter - m_B);
  m_B = m_B - 2 * m_A*parameter;
}

void CostGauss::axisSymmetry(){m_B = - m_B;}

void CostGauss::opposition()
{
  m_A = - m_A;
  m_B = - m_B;
  constant = - constant;
}

//####### minus #######////####### minus #######////####### minus #######//
//####### minus #######////####### minus #######////####### minus #######//


CostGauss CostGauss::minus(CostGauss const& mycost)
{
  CostGauss gaussNew = CostGauss();
  gaussNew.m_A = m_A - mycost.getM_A();
  gaussNew.m_B = m_B - mycost.getM_B();
  gaussNew.constant = constant - mycost.getConstant();

  return(gaussNew);
}




//####### isEqual #######////####### isEqual #######////####### isEqual #######//
//####### isEqual #######////####### isEqual #######////####### isEqual #######//

bool CostGauss::isEqual(CostGauss const& mycost) const
{
  bool res = false;
  res = (mycost.getM_A() == m_A)&&(mycost.getM_B() == m_B)&&(mycost.getConstant() == constant);
  return(res);
}


//####### addL1 #######////####### addL1 #######////####### addL1 #######//
//####### addL1 #######////####### addL1 #######////####### addL1 #######//



void CostGauss::addL1(Point const& pt, Robust const& robust, int slope)
{
  double K = robust.getThreshold();
  double a = robust.getSlope();
  m_B = m_B + slope*pt.w * a;
  constant = constant + K*K - a*K - slope*a*pt.y;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void CostGauss::show() const
{
  //std::cout << "#GAUSS# #A# " << m_A << " #B# "<< m_B << " #C# " << constant;
  //std::cout << " #MIN# "<< minimum();
  //std::cout << " #ARGMIN# "<< arg_minimum() << std::endl;
}


void CostGauss::show(Interval const& inter) const
{
  //std::cout << "#GAUSS# #A# " << m_A << " #B# "<< m_B << " #C# " << constant;
  //std::cout << " #MIN# "<< minimum();
  //double argmin = arg_minimum();
  //if(inter.isInside(argmin) == true){std::cout << " #ARGMIN# "<< argmin << std::endl;}
  //else{std::cout << " #ARGMIN# "<< arg_minimum() << std::endl;}
}




