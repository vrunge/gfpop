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

  //###############################################################
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

  //###############################################################
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

  //###############################################################
  if(edge_ctt == "up")
  {
  }

  //###############################################################
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
    while(currentPiece != NULL)
    {
      currentPiece -> addCoeff(Cost(cost_coeff(pt)), penalty);
      move();
    }
  }

  ///////////////////// CASE K != INF /////////////////////
  ///////////////////// CASE K != INF /////////////////////
  if(K != INFINITY)
  {
    ///DANGER
    double* coeff = new double[3];
    coeff[0] = 0;
    coeff[1] = -a;
    coeff[2] = K;
    Cost slopeLeftCost = Cost(coeff);
    coeff[1] = a;
    Cost slopeRightCost = Cost(coeff);

    ///Interval
    Cost costInter = Cost(cost_coeff(pt));
    Interval new_interval = cost_intervalInterRoots(costInter,K);


    /// Putting the bounds in variables A and B
    double AK = new_interval.geta();
    double BK = new_interval.getb();


    while(currentPiece != NULL)
    {
      currentPiece -> addCoeff(Cost(cost_coeff(pt)), penalty);
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






