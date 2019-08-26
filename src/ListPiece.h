#ifndef LISTPIECE_H
#define LISTPIECE_H

#include "Piece.h"
#include <math.h>

#include <vector>

class ListPiece
{
private:
  Piece* head;
  Piece* currentPosition;
  Piece* lastActivePosition;
  Piece* tail;
  unsigned int lengthList;

public:
  ListPiece();
  ~ListPiece();
  void move();
  void addPiece(Piece* newP);
  void deleteNxtPointAndMove();
  void initializeCurrentPosition();
  unsigned int getLength();

};

#endif // LISTPIECE_H
