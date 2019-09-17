#include "Piece.h"

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


//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//


Track Piece::getTrack()const {return(m_info);}
Interval Piece::getInterval() const {return(m_interval);}
Cost Piece::getCost() const {return(m_cost);}

Cost& Piece::getRefCost(){return(m_cost);}


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




//####### addPointAndPenalty #######////####### addPointAndPenalty #######////####### addPointAndPenalty #######//
//####### addPointAndPenalty #######////####### addPointAndPenalty #######////####### addPointAndPenalty #######//

void Piece::addPointAndPenalty(Point const& pt, double penalty)
{
  //ADD add pt
  addmyConstant(m_cost, penalty);
}



