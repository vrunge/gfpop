// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef LISTPIECE_H
#define LISTPIECE_H

#include "Piece.h"
#include <math.h>

class ListPiece ///DANGER: THERE IS "NO" EMPTY LIST. ALWAYS ADD SOME ELEMENTS
{
public:
  ListPiece(); //create an empty list
  ~ListPiece();

  void addPiece(Piece* P); //at the beginning after (0,0)
  void deleteNxtPiece();

  void initializeCurrentPosition(); //currentPosition = firstPiece
  bool move(); //DANGER we can move only in a non-empty list
  bool isEmpty(); //DANGER we can move only in a non-empty list

private:
  Piece* firstPiece;
  Piece* currentPiece;

};

#endif // LISTPIECE_H
