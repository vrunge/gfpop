#include "ListPiece.h"
#include "Piece.h"
#include <iostream>
#include "stdlib.h"


ListPiece::ListPiece()
{
  lengthList = 0;
  head = new Piece();
  head -> nxt = NULL;
  lastActivePosition = head;
  tail = head;
  currentPosition = head;
}

void ListPiece::addPiece(Piece* newP)
{
  if (lastActivePosition -> nxt == NULL)
  {
    tail = newP;
    lastActivePosition -> nxt = tail;
    lastActivePosition = lastActivePosition->nxt;
  }
  else
  {
    lastActivePosition = lastActivePosition->nxt;
    lastActivePosition = newP;
  }
  lengthList ++;
}

void ListPiece::move()
{
  currentPosition = currentPosition->nxt;
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

void ListPiece::initializeCurrentPosition()
{
  currentPosition = head;
}

void ListPiece::deleteNxtPointAndMove(){
  tail -> nxt = currentPosition -> nxt;
  currentPosition -> nxt = currentPosition->nxt->nxt;
  tail->nxt->nxt = NULL;
  tail = tail->nxt;
  lengthList--;
}



