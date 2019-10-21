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

//##### addNewCurrentPiece #####//////##### addNewCurrentPiece #####//////##### addNewCurrentPiece #####///
//##### addNewCurrentPiece #####//////##### addNewCurrentPiece #####//////##### addNewCurrentPiece #####///

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



// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT ///
// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT ///
// WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY //
// WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY //

// shift_right // shift_right // shift_right // shift_right // shift_right // shift_right //
// shift_right // shift_right // shift_right // shift_right // shift_right // shift_right //

void ListPiece::shift_right(double parameter)
{
}

void ListPiece::shift_left(double parameter)
{
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
  int parentStateLabel = edge.getState1(); ///parentStateLabel = state to associate

  //################
  if(edge_ctt == "null") /// Simple copy of LP_state
  {
    Piece* current = LP_state.head;
    while(current != NULL)
    {
      addNewLastPiece(current -> copy());
      current = current -> nxt;
      ///DANGER : change Track?
    }
  }

  //################
  if(edge_ctt == "std")
  {
    ///find the minimum
    double mini = INFINITY;
    double getmin;
    Piece* current = LP_state.head;
    while(current != NULL)
    {
      getmin = cost_minInterval(current -> getCost(), current -> getInterval());
      if(getmin < mini){mini = getmin;}
      current = LP_state.currentPiece -> nxt;
    }

    ///add one Piece to LP_edges
    ///DANGER : find interval + add the constant = mini
    Piece* onePiece = new Piece();
    addNewLastPiece(onePiece);

  }

  //################
  if(edge_ctt == "up")
  {
  }

  //################
  if(edge_ctt == "down")
  {
    reverse(); ///look at pieces from right to left


    reverse(); ///look at pieces from left to right
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
      currentPiece -> addCoeff(costPt, penalty);
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
      tmpA = currentPiece -> getInterval().geta();
      tmpB = currentPiece -> getInterval().getb();

      if(tmpB <= AK){cas = 0;}
      if(BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;}///priority to AK over BK between tempA and tempB.

      switch(cas)
      {
      case 0 : currentPiece -> addCoeff(slopeLeftCost, penalty);
        break;
      case 1 : currentPiece -> addCoeff(slopeRightCost, penalty);
        break;
      case 2 : currentPiece -> addCoeff(costPt, penalty);
        break;
      case 3 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = currentPiece -> copy();
        addCurrentPiecePlus1(new_piece);
        ///adding costPt on the left
        currentPiece -> addCoeff(costPt, penalty);
        ///changing interval bounds
        currentPiece -> setIntervalB(BK);
        new_piece -> setIntervalA(BK);
        break;
      }

      case 4 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = new Piece(currentPiece);
        addCurrentPiecePlus1(new_piece);
        ///adding slopeLeftCost on the left
        currentPiece -> addCoeff(slopeLeftCost, penalty);
        ///changing interval bounds
        currentPiece -> setIntervalB(AK);
        new_piece -> setIntervalA(AK);
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


/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::setUniquePieceCostToInfinity()
{
  head -> getRefCost().constant = INFINITY;
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






