#include "ListPiece.h"
#include "Piece.h"
#include <iostream>
#include "stdlib.h"


ListPiece::ListPiece()
{
  head = NULL;
  currentPiece = NULL;
  lastPiece = NULL;
}

ListPiece::~ListPiece()
{
  while(head != NULL)
  {
    Piece* pieceToDelete = head;
    head = head -> nxt;
    delete(pieceToDelete);
  }
}

//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///
//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///

void ListPiece::setUniquePieceCostToInfinity()
{
  head -> m_cost.constant = INFINITY;
}

//##### reset #####//////##### reset #####//////##### reset #####///
//##### reset #####//////##### reset #####//////##### reset #####///

void ListPiece::reset()
{
  while(head != NULL)
  {
    Piece* pieceToDelete = head;
    head = head -> nxt;
    delete(pieceToDelete);
  }
  currentPiece = NULL;
  lastPiece = NULL;
}

//##### reverse #####//////##### reverse #####//////##### reverse #####///
//##### reverse #####//////##### reverse #####//////##### reverse #####///

void ListPiece::reverse()
{
  lastPiece = head;

  Piece* prev =  NULL;
  Piece* current = head;
  Piece* next = current;

  while(current != NULL)
  {
    next  = current -> nxt;
    current -> nxt = prev; /// new nxt
    prev = current;
    current = next;
  }

  head = prev;
  currentPiece = head;
}

//##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####///
//##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####///

void ListPiece::addCurrentPiecePlus1(Piece* newPiece)
{
  newPiece -> nxt = currentPiece -> nxt;
  currentPiece -> nxt = newPiece;
}

//##### addNewLastPiece #####//////##### addNewLastPiece #####//////##### addNewLastPiece #####///
//##### addNewLastPiece #####//////##### addNewLastPiece #####//////##### addNewLastPiece #####///

void ListPiece::addNewLastPiece(Piece* newPiece)
{
  lastPiece -> nxt = newPiece;
  lastPiece = newPiece;
}

//##### move #####//////##### move #####//////##### move #####///
//##### move #####//////##### move #####//////##### move #####///

void ListPiece::move()
{
  currentPiece = currentPiece -> nxt;
}

//##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####///
//##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####///

void ListPiece::initializeCurrentPiece()
{
  currentPiece = head;
}

//##### copy #####//////##### copy #####//////##### copy #####///
//##### copy #####//////##### copy #####//////##### copy #####///

void ListPiece::copy(ListPiece const& LP_edge)
{
  Piece* current = LP_edge.head;
  head = current;

  while(current != NULL)
  {
    currentPiece = current -> copy(); //copy content in Piece
    currentPiece -> nxt = current -> nxt; //copy nxt pointer
    current = current -> nxt;
  }
  lastPiece = currentPiece;
}


//##### shift #####//////##### shift #####//////##### shift #####///
//##### shift #####//////##### shift #####//////##### shift #####///

void ListPiece::shift(double parameter)
{
  currentPiece = head;
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    Interval inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interShift(inter.geta(), parameter));
    currentPiece -> m_interval.setb(cost_interShift(inter.getb(), parameter));
    ///MOVE Cost
    cost_shift(currentPiece -> m_cost, parameter);
    currentPiece = currentPiece -> nxt;
  }
}


//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///
//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///

void ListPiece::expDecay(double gamma)
{
  currentPiece = head;
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    Interval inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interExpDecay(inter.geta(), gamma));
    currentPiece -> m_interval.setb(cost_interExpDecay(inter.getb(), gamma));
    ///MOVE Cost
    cost_expDecay(currentPiece -> m_cost, gamma);
    currentPiece = currentPiece -> nxt;
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///
//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///

void ListPiece::LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel)
{
  reset(); /// build a new LP_edges from scratch

  /// 4 types of edges : null, std, up, down
  ///
  /// EDGE PARAMETERS
  ///
  std::string edge_ctt = edge.getConstraint();
  double edge_parameter = edge.getParameter();
  double edge_beta = edge.getBeta();
  unsigned int parentState = edge.getState1(); ///parentState = state to associate

  //################
  if(edge_ctt == "null") /// Simple copy of LP_state
  {
    copy(LP_state);
    if(edge_parameter < 1){expDecay(edge_parameter);}
  }

  //################
  if(edge_ctt == "std")
  {
    ///find the minimum
    double mini = INFINITY;
    double getmin;
    unsigned int counter = 0;
    unsigned int myCounter;

    Piece* current = LP_state.head;
    while(current != NULL)
    {
      counter = counter + 1;
      getmin = cost_minInterval(current -> m_cost, current -> m_interval);
      if(getmin < mini){mini = getmin; myCounter = counter;}
      current = current -> nxt;
    }

    ///add onePiece to LP_edges
    Piece* onePiece = new Piece();
    onePiece ->addCostAndPenalty(Cost(), mini + edge_beta); /// Cost() = 0
    onePiece -> m_info = Track(newLabel, parentState, myCounter);
    addNewLastPiece(onePiece);

  }

  //################
  if(edge_ctt == "up")
  {
    operator_up(LP_state, newLabel, parentState, true);
    if(edge_parameter < 1){shift(edge_parameter);}
  }

  //################
  if(edge_ctt == "down")
  {
    reverse(); ///look at pieces from right to left
    operator_up(LP_state, newLabel, parentState, false);
    reverse(); ///look at pieces from left to right
    if(edge_parameter < 1){shift(-edge_parameter);}
  }

}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void ListPiece::LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt)
{
  /// get edge data ///
  double K = edge.getKK();
  double a = edge.getAA();
  double penalty = edge.getBeta();

  initializeCurrentPiece();

  ///////////////////// CASE K == INF /////////////////////
  ///////////////////// CASE K == INF /////////////////////
  if(K == INFINITY)
  {
    Cost costPt = Cost(cost_coeff(pt));
    while(currentPiece != NULL)
    {
      currentPiece -> addCostAndPenalty(costPt, penalty);
      move();
    }
  }

  ///////////////////// CASE K != INF /////////////////////
  ///////////////////// CASE K != INF /////////////////////
  if(K != INFINITY)
  {
    /// INTIALIZATION
    ///Interval
    Cost costInter = Cost(cost_coeff(pt));
    Interval new_interval = cost_intervalInterRoots(costInter, K);

    /// Putting the bounds in variables AK and BK
    double AK = new_interval.geta();
    double BK = new_interval.getb();

    double* coeff = new double[3];
    coeff[0] = 0;
    coeff[1] = -a;
    coeff[2] = K + a * AK;
    Cost slopeLeftCost = Cost(coeff);  /// LEFT y = -ax + K + a * AK
    coeff[1] = a;
    coeff[2] = K  - a * BK;
    Cost slopeRightCost = Cost(coeff);  /// RIGHT y = ax + K - a * BK
    Cost costPt = Cost(cost_coeff(pt));  /// CENTER pt

    /// bounds of the tmp Piece
    double tmpA = 0;
    double tmpB = 0;

    int cas = 0;

    while(currentPiece != NULL)
    {
      tmpA = currentPiece -> m_interval.geta();
      tmpB = currentPiece -> m_interval.getb();

      if(tmpB <= AK){cas = 0;}
      if(BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;}///priority to AK over BK between tempA and tempB.

      switch(cas)
      {
      case 0 : currentPiece -> addCostAndPenalty(slopeLeftCost, penalty);
        break;
      case 1 : currentPiece -> addCostAndPenalty(slopeRightCost, penalty);
        break;
      case 2 : currentPiece -> addCostAndPenalty(costPt, penalty);
        break;
      case 3 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = currentPiece -> copy();
        addCurrentPiecePlus1(new_piece);
        ///adding costPt on the left
        currentPiece -> addCostAndPenalty(costPt, penalty);
        ///changing interval bounds
        currentPiece -> m_interval.setb(BK);
        new_piece -> m_interval.seta(BK);
        break;
      }

      case 4 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = new Piece(currentPiece);
        addCurrentPiecePlus1(new_piece);
        ///adding slopeLeftCost on the left
        currentPiece -> addCostAndPenalty(slopeLeftCost, penalty);
        ///changing interval bounds
        currentPiece -> m_interval.setb(AK);
        new_piece -> m_interval.seta(AK);
        break;
      }
      }
    move();
    }
  }


}



//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///
//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///

void ListPiece::LP_ts_Minimization(ListPiece const& LP_edge)
{


}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### operator_up #####//////##### operator_up #####//////##### operator_up #####///
//##### operator_up #####//////##### operator_up #####//////##### operator_up #####///

void ListPiece::operator_up(ListPiece const& LP_edge, unsigned int newLabel, unsigned int parentState, bool upDirection)
{
  Piece* tmp = LP_edge.head;
  double currentValue;
  head = new Piece();
  unsigned int counter = 1;

  Track trackUp = Track(newLabel, parentState, counter);
  head -> m_info.setTrack(trackUp); ///set Track


  while(tmp != NULL)
  {
    currentPiece -> paste(tmp, currentValue);
    tmp = tmp -> nxt;
    counter = counter + 1;
  }

}


/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::show()
{
  initializeCurrentPiece();

  while(currentPiece != NULL)
  {
    currentPiece -> show();
    move();
  }
}

