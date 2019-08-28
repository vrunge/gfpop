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


void ListPiece::addPointAndPenalty(Point const& pt, Edge const& edge)
{
  double K = edge.getKK();
  double a = edge.getAA();
  double penalty = edge.getBeta();

  initializeCurrentPiece();

  while(currentPiece != NULL)
  {
    currentPiece -> addPointAndPenalty(pt, penalty);
    move();
  }

}

/////////////////////////////////////////


ListPiece ListPiece::edgeConstraintLP(Edge const& edge, int newLabel, Bound const& bound)
{
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






