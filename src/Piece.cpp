/* -*- compile-command: "R CMD INSTALL .." -*- */
#include "Piece.h"
#include"ExternFunctions.h"

#include <stdlib.h>
#include<iostream>

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

Piece* Piece::copy(){return(new Piece(this));}

//####### addCostAndPenalty #######////####### addCostAndPenalty #######////####### addCostAndPenalty #######//
//####### addCostAndPenalty #######////####### addCostAndPenalty #######////####### addCostAndPenalty #######//

void Piece::addCostAndPenalty(Cost const& cost, double penalty)
{
  m_cost.m_A = m_cost.m_A + cost.m_A;
  m_cost.m_B = m_cost.m_B + cost.m_B;
  m_cost.constant = m_cost.constant + cost.constant + penalty;
}


//####### intervalMinLessUp #######////####### intervalMinLessUp #######////####### intervalMinLessUp #######//
//####### intervalMinLessUp #######////####### intervalMinLessUp #######////####### intervalMinLessUp #######//


Interval Piece::intervalMinLessUp(double bound, double currentValue, bool constPiece)
{
  //
  // (bound, currentValue) VERSUS (argmin, min)
  //

  Interval response = Interval(); /// Interval = (INFINITY, INFINITY)
  double mini = cost_minInterval(m_cost, m_interval); /// It was: double mini = cost_min(m_cost);

  if(currentValue > mini) /// otherwise currentValue constant doesn't intersect Piece cost
  {
    double argmini = cost_argmin(m_cost);
    if(bound < argmini) /// otherwise currentValue constant doesn't intersect Piece cost
    {
      if(constPiece == true)
      {
        double* coeff = new double[3];
        coeff[0] = m_cost.m_A;
        coeff[1] = m_cost.m_B;
        coeff[2] = m_cost.constant;
        Cost costInter = Cost(coeff);
        response.seta(cost_intervalInterRoots(costInter, currentValue).geta());
        delete[] coeff;
      }
      else
      {
        response.seta(bound); /// i.e. cost_eval(bound) == mini : continuity condition
      }
    response.setb(argmini);
    }
  }
  else if(currentValue == mini)
  {
    response.seta(bound); /// i.e. cost_eval(bound) == mini : continuity condition
    response.setb(m_interval.getb());
  }
  return(response);
}



//####### intervalMinLessDw #######////####### intervalMinLessDw #######////####### intervalMinLessDw #######//
//####### intervalMinLessDw #######////####### intervalMinLessDw #######////####### intervalMinLessDw #######//

Interval Piece::intervalMinLessDw(double bound, double currentValue, bool constPiece)
{
  //
  // (bound, currentValue) VERSUS (argmin, min)
  //

  Interval response = Interval(); /// Interval = (INFINITY, INFINITY)
  double mini = cost_minInterval(m_cost, m_interval); /// It was: double mini = cost_min(m_cost);

  if(currentValue > mini) /// otherwise currentValue constant doesn't intersect Piece cost
  {
    double argmini = cost_argmin(m_cost);
    if(bound > argmini) /// otherwise currentValue constant doesn't intersect Piece cost
    {
      if(constPiece == true)
      {
        double* coeff = new double[3];
        coeff[0] = m_cost.m_A;
        coeff[1] = m_cost.m_B;
        coeff[2] = m_cost.constant;
        Cost costInter = Cost(coeff);
        response.setb(cost_intervalInterRoots(costInter, currentValue).getb());
        delete[] coeff;
      }
      else
      {
        response.setb(bound); /// i.e. cost_eval(bound) == mini : continuity condition
      }
    response.seta(argmini);
    }
  }
  else if(currentValue == mini)
  {
    response.seta(m_interval.geta());
    response.setb(bound); /// i.e. cost_eval(bound) == mini : continuity condition
  }
  return(response);
}

//####### pastePieceUp #######// //####### pastePieceUp #######// //####### pastePieceUp #######//
//####### pastePieceUp #######// //####### pastePieceUp #######// //####### pastePieceUp #######//


Piece* Piece::pastePieceUp(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack)
{
  Piece* BUILD = this;

  /// decreasingInterval = (a,b)
  /// Q -> m_interval = (m_a,m_b)

  if(decrInter.isEmpty())
  {
    BUILD -> m_interval.setb(NXTPiece -> m_interval.getb());
  }
  else
  {
    BUILD -> m_interval.setb(decrInter.geta());    ///if a > m_a, no change otherwise

    ///ADD of the PIECE NXTPiece (troncated)
    if(BUILD -> m_interval.isEmpty()) ///if BUILD empty
    {
      BUILD -> m_interval.setb(decrInter.getb());
      BUILD -> m_cost = NXTPiece -> m_cost;
      BUILD -> m_info.setTrack(newTrack);
    }
    else
    {
      Piece* NewQ = new Piece(newTrack, decrInter, NXTPiece -> m_cost);
      BUILD -> nxt = NewQ;
      BUILD = NewQ;
    }

    if(!((NXTPiece -> nxt == NULL) && (decrInter.getb() == NXTPiece -> m_interval.getb())))
    {
      double outputValue = cost_eval(NXTPiece -> m_cost, decrInter.getb());
      Piece* PieceOut = new Piece(newTrack, Interval(decrInter.getb(), NXTPiece -> m_interval.getb()), Cost());
      addConstant(PieceOut -> m_cost, outputValue);
      BUILD -> nxt = PieceOut;
      BUILD = PieceOut;
    }
  }

  return(BUILD);
}



//####### pastePieceDw #######// //####### pastePieceDw #######// //####### pastePieceDw #######//
//####### pastePieceDw #######// //####### pastePieceDw #######// //####### pastePieceDw #######//


Piece* Piece::pastePieceDw(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack)
{
  Piece* BUILD = this;

  /// decreasingInterval = (a,b)
  /// NXTPiece -> m_interval = (m_a,m_b)

  if(decrInter.isEmpty())
  {
    BUILD -> m_interval.seta(NXTPiece -> m_interval.geta());
  }
  else
  {
    BUILD -> m_interval.seta(decrInter.getb());    ///if a > m_a, no change otherwise

    ///ADD of the PIECE NXTPiece (troncated)
    if(BUILD -> m_interval.isEmpty()) ///if BUILD empty
    {
      BUILD -> m_interval.seta(decrInter.geta());
      BUILD -> m_cost = NXTPiece -> m_cost;
      BUILD -> m_info.setTrack(newTrack);
    }
    else
    {
      Piece* NewQ = new Piece(newTrack, decrInter, NXTPiece -> m_cost);
      BUILD -> nxt = NewQ;
      BUILD = NewQ;
    }

    if(!((NXTPiece -> nxt == NULL) && (decrInter.geta() == NXTPiece -> m_interval.geta())))
    {
      double outputValue = cost_eval(NXTPiece -> m_cost, decrInter.geta());
      Piece* PieceOut = new Piece(newTrack, Interval(NXTPiece -> m_interval.geta(), decrInter.geta()), Cost());
      addConstant(PieceOut -> m_cost, outputValue);
      BUILD -> nxt = PieceOut;
      BUILD = PieceOut;
    }
  }
  return(BUILD);
}


//####### pieceGenerator #######// //####### pieceGenerator #######// //####### pieceGenerator #######//
//####### pieceGenerator #######// //####### pieceGenerator #######// //####### pieceGenerator #######//

Piece* Piece::pieceGenerator(Piece* Q1, Piece* Q2, int Bound_Q2_Minus_Q1, double M)
{
  Piece* BUILD = this; // = Q12
  double zero = 0;

  // INFORMATION interToPaste
  // Interval interToPaste = interval on which we build BUILD = right BUILD TO min(right Q1, right Q2)
  Interval interToPaste = Interval();
  interToPaste.seta(BUILD -> m_interval.getb()); ///enter Piece BUILD
  if(Bound_Q2_Minus_Q1 == -1)
         {interToPaste.setb(Q2 -> m_interval.getb());
    }else{interToPaste.setb(Q1 -> m_interval.getb());}

  // INFORMATION interRoots
  // Interval interRoots (Q1 - Q2)
  Cost costDiff = minusCost(Q1 -> m_cost, Q2 -> m_cost);
  Interval interRoots = cost_intervalInterRoots(costDiff, zero);
  //std::cout << "interRoots ";
  //interRoots.show();
  //showCost(costDiff);

  //// INFORMATION change
  // int change = 0, 1 or 2 change-points
  // if bounds of interRoots close to bounds of interToPaste > 1e-12
  unsigned int change = 0;
  if((interRoots.geta() > interToPaste.geta() + 1e-12)&&(interRoots.geta() + 1e-12 < interToPaste.getb())){change = change + 1;}
  if((interRoots.getb() > interToPaste.geta() + 1e-12)&&(interRoots.getb() + 1e-12 < interToPaste.getb())){change = change + 1;}

  ///Security steps: length interRoots very small < 1e-12
  if(interRoots.getb() - interRoots.geta() < 1e-12)
  {
    change = 0;
    interRoots.seta(interToPaste.geta());
    interRoots.setb(interToPaste.getb());
  }

  // CONSTRUCTION
  int Q2_Minus_Q1;  ///Sign of Q2 - Q1

  switch(change)
  {
    case 0 : /// IF WE ADD 0 PIECE
    {
      BUILD = BUILD -> piece0(Q1, Q2, interToPaste, Q2_Minus_Q1);
      break;
    }
    case 1 : /// IF WE ADD 1 PIECE
    {
      BUILD = BUILD -> piece1(Q1, Q2, interToPaste, interRoots, Q2_Minus_Q1);
      break;
    }
    case 2 : /// IF WE ADD 2 PIECES
    {
      BUILD = BUILD -> piece2(Q1, Q2, interToPaste, interRoots, Q2_Minus_Q1);
      break;
    }
  }

  //CONSTRUCTION outPiece
  ///Need of a last "OUT" Piece if we have reached the end of the Piece Q1 or Piece Q2
  if((((Q2_Minus_Q1 == 1) && (Bound_Q2_Minus_Q1 >= 0)) || ((Q2_Minus_Q1 == -1) && (Bound_Q2_Minus_Q1 <= 0)))
       && (interToPaste.getb() != M))
  {
    Piece* outPiece = new Piece();
    outPiece -> m_interval = Interval(interToPaste.getb(), interToPaste.getb());
    BUILD -> nxt = outPiece;
    BUILD = outPiece;
  }
  return(BUILD);
}


//####### piece0 #######// //####### piece0 #######// //####### piece0 #######//
//####### piece0 #######// //####### piece0 #######// //####### piece0 #######//

Piece* Piece::piece0(Piece* Q1, Piece* Q2, Interval interToPaste, int& Q2_Minus_Q1)
{
  Piece* BUILD = this;

  /// Possible inversion => test Q2_Minus_Q1 at centerPoint
  double centerPoint = interToPaste.internPoint();
  Cost costDiff = minusCost(Q2 -> m_cost, Q1 -> m_cost);
  Q2_Minus_Q1 = signValue(cost_eval(costDiff, centerPoint));
  bool test = false;

  if (BUILD -> m_interval.isEmpty() == true) /// IF BUILD interval = empty
  {
    BUILD -> m_interval.setb(interToPaste.getb());
    if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
    if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}
  }
  else
  {
    ///SECURITY step
    Cost testCost = BUILD -> m_cost;
    if(Q2_Minus_Q1 == 1){test = isEqual(testCost, Q1 -> m_cost);}
    if(Q2_Minus_Q1 == -1){test = isEqual(testCost, Q2 -> m_cost);}

    if (test == true) ///Prolongation
    {
      BUILD -> m_interval.setb(interToPaste.getb());
      if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}
    }
    else ///pb with the cost -> we stop BUILD interval at interToPaste left -> we create a new piece
    {
      //CONSTRUCTION newPiece
      BUILD -> m_interval.setb(interToPaste.geta());
      Piece* newPiece = new Piece();
      newPiece -> m_interval = interToPaste;
      if(Q2_Minus_Q1 == 1){newPiece -> m_cost = Q1 -> m_cost; newPiece -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){newPiece -> m_cost = Q2 -> m_cost; newPiece -> m_info = Q2 -> m_info;}
      BUILD -> nxt = newPiece;
      BUILD = newPiece;
    }
  }
  return(BUILD);
}

//####### piece1 #######// //####### piece1 #######// //####### piece1 #######//
//####### piece1 #######// //####### piece1 #######// //####### piece1 #######//

Piece* Piece::piece1(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1)
{
  Piece* BUILD = this;
  //PROLONGATION theChangePoint
  double theChangePoint;
  if(interToPaste.geta() < interRoots.geta()){theChangePoint = interRoots.geta();}
                                         else{theChangePoint = interRoots.getb();}

  //// FIND the winner on the new piece
  // centerPoint = centre (left interToPaste, right theChangePoint)
  double centerPoint = Interval(interToPaste.geta(), theChangePoint).internPoint();
  Cost costDiff = minusCost(Q2 -> m_cost, Q1 -> m_cost);
  Q2_Minus_Q1 = signValue(cost_eval(costDiff, centerPoint));

  if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
  if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}

  BUILD -> m_interval.setb(theChangePoint);

  //CONSTRUCTION newPiece
  Piece* newPiece = new Piece();
  newPiece -> m_interval = Interval(theChangePoint, interToPaste.getb());

  centerPoint = newPiece -> m_interval.internPoint();
  costDiff = minusCost(Q2 -> m_cost, Q1 -> m_cost);
  Q2_Minus_Q1 = signValue(cost_eval(costDiff, centerPoint));

  if(Q2_Minus_Q1 == 1){newPiece -> m_cost = Q1 -> m_cost; newPiece -> m_info = Q1 -> m_info;}
  if(Q2_Minus_Q1 == -1){newPiece -> m_cost = Q2 -> m_cost; newPiece -> m_info = Q2 -> m_info;}
  BUILD -> nxt = newPiece;
  BUILD = newPiece;

  return(BUILD);
}




//####### piece2 #######// //####### piece2 #######// //####### piece2 #######//
//####### piece2 #######// //####### piece2 #######// //####### piece2 #######//

Piece* Piece::piece2(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1)
{
  Piece* BUILD = this;

  //PROLONGATION theChangePoint
  //FIND the winner on the newpiece1 (the central piece defined on interval interRoots)
  double centerPoint = interRoots.internPoint();
  Cost costDiff = minusCost(Q2 -> m_cost, Q1 -> m_cost);
  Q2_Minus_Q1 = signValue(cost_eval(costDiff, centerPoint));///INVERSION!!!

  Q2_Minus_Q1 = -Q2_Minus_Q1; ///INVERSION!!!
  if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
  if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}
  BUILD -> m_interval.setb(interRoots.geta());

  //CONSTRUCTION newPiece1
  Q2_Minus_Q1 = -Q2_Minus_Q1; ///INVERSION!!!
  Piece* newPiece1 = new Piece();
  newPiece1 -> m_interval = interRoots;
  if(Q2_Minus_Q1 == 1){newPiece1 -> m_cost = Q1 -> m_cost; newPiece1 -> m_info = Q1 -> m_info;}
  if(Q2_Minus_Q1 == -1){newPiece1 -> m_cost = Q2 -> m_cost; newPiece1 -> m_info = Q2 -> m_info;}
  BUILD -> nxt = newPiece1;
  BUILD = newPiece1;

  Q2_Minus_Q1 = -Q2_Minus_Q1;  ///INVERSION!!!
  //CONSTRUCTION newPiece2
  Piece* newPiece2 = new Piece();
  newPiece2 -> m_interval = Interval(interRoots.getb(), interToPaste.getb());
  if(Q2_Minus_Q1 == 1){newPiece2 -> m_cost = Q1 -> m_cost; newPiece2 -> m_info = Q1 -> m_info;}
  if(Q2_Minus_Q1 == -1){newPiece2 -> m_cost = Q2 -> m_cost; newPiece2 -> m_info = Q2 -> m_info;}
  BUILD -> nxt = newPiece2;
  BUILD = newPiece2;

  return(BUILD);
}


//####### get_min_argmin_label_state_position #######// //####### get_min_argmin_label_state_position #######// //####### get_min_argmin_label_state_position #######//
//####### get_min_argmin_label_state_position #######// //####### get_min_argmin_label_state_position #######// //####### get_min_argmin_label_state_position #######//

void Piece::get_min_argmin_label_state_position(double* response)
{
  response[0] = cost_minInterval(this -> m_cost, this -> m_interval);
  response[1] = cost_argminBacktrack(this -> m_cost, this -> m_interval);
  response[2] = this -> m_info.getLabel();
  response[3] = this -> m_info.getState();
  response[4] = this -> m_info.getPosition();
}



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


void Piece::show()
{
  Piece* tmp = this;
  //if(tmp == NULL){std::cout << "#NULL EMPTY POINTER# "<< std::endl;}
  //else
  {
    //std::cout << "          " << tmp;
    //std::cout << " #LABEL# "<< tmp -> m_info.getLabel() << " #STATE# " <<  tmp -> m_info.getState() << " POSITION " << tmp -> m_info.getPosition() << " ";
    //std::cout << " #INTERVAL# "<< tmp -> m_interval.geta() << " to " << tmp -> m_interval.getb() << " ";
    //showCost(tmp -> m_cost);
  }
  if(tmp -> nxt != NULL){tmp -> nxt -> show();}
}



