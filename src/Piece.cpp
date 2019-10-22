#include "Piece.h"
#include"ExternFunctions.h"

#include <math.h>
#include <stdlib.h>
#include<iostream>
#include<typeinfo>

#include <fstream> ///write in a file
#include<vector>

Piece::Piece(){m_info = Track(); m_interval = Interval(); m_cost = Cost(); nxt = NULL;}
Piece::Piece(Track const& info, Interval const& inter, Cost const& cost)
{
  m_info = info;
  m_interval = inter;
  m_cost = cost;
  nxt = NULL;
}

Piece::Piece(const Piece* piece)
{
  m_info = piece -> m_info;
  m_interval = piece -> m_interval;
  m_cost = piece -> m_cost;
  nxt = NULL;
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Piece::~Piece()
{
  delete(nxt);
  nxt = NULL;
}



//####### copy #######////####### copy #######////####### copy #######//
//####### copy #######////####### copy #######////####### copy #######//

Piece* Piece::copy()
{
  return(new Piece(this));
}



//####### getMin #######////####### getMin #######////####### getMin #######//
//####### getMin #######////####### getMin #######////####### getMin #######//

double Piece::getMin()
{
  return(cost_minInterval(m_cost, m_interval));
}



//####### addCostAndPenalty #######////####### addCostAndPenalty #######////####### addCostAndPenalty #######//
//####### addCostAndPenalty #######////####### addCostAndPenalty #######////####### addCostAndPenalty #######//

void Piece::addCostAndPenalty(Cost const& cost, double penalty)
{
  m_cost.m_A = m_cost.m_A + cost.m_A;
  m_cost.m_B = m_cost.m_B + cost.m_B;
  m_cost.constant = m_cost.constant + cost.constant + penalty;
}


//####### paste #######////####### paste #######////####### paste #######//
//####### paste #######////####### paste #######////####### paste #######//

void Piece::paste(Piece* tmp, double currentValue)
{
  ///INFORMATION
  double argmini = cost_argmin(tmp -> m_cost);
  double mini = cost_minInterval(tmp -> m_cost, tmp -> m_interval);
  Interval inter = cost_intervalInterRoots(tmp -> m_cost, currentValue);



}



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


void Piece::show()
{
  Piece* tmp = this;
  if(tmp == NULL){std::cout << "#NULL EMPTY POINTER# "<< std::endl;}
  else
  {
    std::cout << tmp;
    std::cout << " #LABEL# "<< tmp -> m_info.getLabel() << " #STATE# " <<  tmp -> m_info.getState() << " POSITION " << tmp -> m_info.getPosition() << " ";
    std::cout << "#INTERVAL# "<< tmp -> m_interval.geta() << " -- " << tmp -> m_interval.getb() << " ";
    showCost(tmp ->m_cost);
    std::cout << std::endl;
  }
}



