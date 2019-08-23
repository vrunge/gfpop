// MIT License
// Copyright (c) 2019 Vincent Runge

#include "ListPiece.h"
#include <iostream>
#include "math.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

ListPiece::ListPiece()
{
  firstPiece = new Piece();
  currentPiece = firstPiece;
  //the first Piece = an empty Piece : we consider Piece->nxt in functions
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

ListPiece::~ListPiece()
{
  while(firstPiece != NULL)
  {
    Piece* PieceToDelete = firstPiece;
    firstPiece = firstPiece -> nxt;
    delete(PieceToDelete);
  }
}

//####### addPiece #######////####### addPiece #######////####### addPiece #######//
//####### addPiece #######////####### addPiece #######////####### addPiece #######//

void ListPiece::addPiece(Piece* P)
{
  P -> nxt = firstPiece -> nxt;
  firstPiece -> nxt = P;
}


//####### deleteNxtPoint #######////####### deleteNxtPoint #######////####### deleteNxtPoint #######//
//####### deleteNxtPoint #######////####### deleteNxtPoint #######////####### deleteNxtPoint #######//


void ListPiece::deleteNxtPiece()
{
  if(currentPiece -> nxt != NULL)
  {
    Piece* PieceToDelete = currentPiece -> nxt;
    currentPiece -> nxt = currentPiece -> nxt -> nxt;
    delete(PieceToDelete);
  }
}



//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//
//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//

void ListPiece::initializeCurrentPosition()
{currentPiece = firstPiece;}


//####### move #######////####### move #######////####### move #######//
//####### move #######////####### move #######////####### move #######//

bool ListPiece::move()
{
  bool res = true;
  if(currentPiece -> nxt != NULL && currentPiece -> nxt -> nxt != NULL){currentPiece = currentPiece -> nxt;}
    else{res = false;}

  return(res);
}

//####### isEmpty #######////####### isEmpty #######////####### isEmpty #######//
//####### isEmpty #######////####### isEmpty #######////####### isEmpty #######//

bool ListPiece::isEmpty()
{
  bool res = false;
  if(firstPiece -> nxt == NULL){res = true;}
  return(res);
}





