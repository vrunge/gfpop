#include "ListPiece.h"
#include "Piece.h"
#include <iostream>
#include "stdlib.h"


ListPiece::ListPiece()
{
  lengthList = 0;
  head = new Piece();
  head -> nxt = NULL;
  lastActivePiece = head;
  tail = head;
  currentPiece = head;
}

void ListPiece::addPiece(Piece* newP)
{
  if (lastActivePiece -> nxt == NULL)
  {
    tail = newP;
    tail -> nxt = NULL;
    lastActivePiece -> nxt = tail;
    lastActivePiece = tail;
  }
  else
  {
    lastActivePiece = lastActivePiece->nxt;
    lastActivePiece = newP;
  }

  lengthList ++;
}

void ListPiece::move()
{
  currentPiece = currentPiece->nxt;
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

void ListPiece::deleteNxtPieceAndMove(){
  tail -> nxt = currentPiece -> nxt;
  currentPiece -> nxt = currentPiece->nxt->nxt;
  tail = tail->nxt;
  tail->nxt = NULL;
  lengthList--;
}

/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::LP_edges_addPointAndPenalty(Point const& pt, Edge const& edge)
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


ListPiece ListPiece::LP_edges_constraint(Edge const& edge, unsigned int t)
{
  ListPiece LP_edgesNew = ListPiece();
  return(LP_edgesNew);
}


void ListPiece::addConstant(double myconst)
{
}


void ListPiece::show()
{
  initializeCurrentPiece();

  while(currentPiece != NULL)
  {
    currentPiece -> show();
    move();
  }
}






