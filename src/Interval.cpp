#include "Interval.h"

#include<iostream>

Interval::Interval() : m_a(INFINITY), m_b(INFINITY){} ///empty interval = (+INF, +INF)
Interval::Interval(double a, double b) : m_a(a), m_b(b){}

bool Interval::isEmpty() const {return(m_a >= m_b);}

Interval Interval::intersection(Interval const& inter) const
{
  Interval response = *this;
  if((inter.m_a < m_b) && (m_a < inter.m_b))
  {
    response.m_a = std::max(m_a,inter.m_a);
    response.m_b = std::min(m_b,inter.m_b);
  }
  else
  {
    response.m_a = INFINITY;
    response.m_b = INFINITY;
  }
  return(response);
}

void Interval::seta(double a){m_a = a;}
void Interval::setb(double b){m_b = b;}
double Interval::geta() const{return(m_a);}
double Interval::getb() const{return(m_b);}

bool Interval::isInside(double x) const {return((m_a <= x) && (x <= m_b));}

double Interval::internPoint() const
{
  double thePoint = INFINITY;
  if((m_a == -INFINITY) && (m_b == INFINITY)){thePoint = 0;}
  if((m_a == -INFINITY) && (m_b < INFINITY)){thePoint = m_b - 1;}
  if((-INFINITY < m_a) && (m_b == INFINITY)){thePoint = m_a + 1;}
  if((-INFINITY < m_a) && (m_b < INFINITY)){thePoint = (m_a + 2*m_b)/3;}
  return(thePoint);
}


void Interval::show() const
{
  //std::cout << "INTERVAL ## " << m_a << " -- " << m_b << std::endl;
}
