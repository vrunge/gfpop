#include "ListPiece.h"
#include "Piece.h"
#include <iostream>
#include "stdlib.h"


ListPiece::ListPiece()
{
  lengthList = 0;
  head = new Piece();
  head -> nxt = NULL;
  currentPiece = head;
}

void ListPiece::addPiece(Piece* newP)
{
  currentPiece = currentPiece -> nxt;
  currentPiece = newP;
  lengthList = lengthList + 1;
}

void ListPiece::move()
{
  currentPiece = currentPiece -> nxt;
}

unsigned int ListPiece::getLength()
{
  return lengthList;
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

void ListPiece::initializeCurrentPiece()
{
  currentPiece = head;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///
//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///

ListPiece ListPiece::LP_edges_constraint(Edge const& edge, unsigned int t)
{
  ListPiece LP_edgesNew = ListPiece();
  initializeCurrentPiece();

  while(currentPiece != NULL)
  {
    move();
  }
  return(LP_edgesNew);
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






