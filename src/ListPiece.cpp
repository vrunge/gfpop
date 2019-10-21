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
  if(edge_ctt == "null")
  {
  }
  //###############################################################
  if(edge_ctt == "std")
  {
  }
  //###############################################################
  if(edge_ctt == "down")
  {
  }

  //###############################################################
  if(edge_ctt == "up")
  {
  }

  while(currentPiece != NULL)
  {
    move();
  }
}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void ListPiece::LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt)
{
  double K = edge.getKK();
  double a = edge.getAA();
  double penalty = edge.getBeta();
  initializeCurrentPiece();

  ///////////////////// CASE K == INF /////////////////////
  if(K == INFINITY)
  {
    while(currentPiece != NULL)
    {
      currentPiece -> addPointAndPenalty(pt, penalty);
      move();
    }
  }

  ///////////////////// CASE K != INF /////////////////////
  if(K != INFINITY)
  {
    double* coeff = new double[3];
    coeff[0] = 0;
    coeff[1] = a;
    coeff[2] = K;
    Cost cost = Cost(coeff);

    while(currentPiece != NULL)
    {
      currentPiece -> addPointAndPenalty(pt, penalty);
      move();
    }
  }


}

/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::addConstant(double myconst)
{

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






