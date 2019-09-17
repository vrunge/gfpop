#ifndef LISTPIECE_H
#define LISTPIECE_H

#include "Piece.h"
#include "Edge.h"
#include "ExternFunctions.h"

#include <math.h>

class ListPiece
{
private:
  Piece* head;
  Piece* currentPiece;
  Piece* lastActivePiece;
  Piece* tail;
  unsigned int lengthList;

public:
  ListPiece();
  ~ListPiece();
  void move();
  void addPiece(Piece* newP);
  void deleteNxtPieceAndMove();
  void initializeCurrentPiece();
  unsigned int getLength();

  ///////  3 OPERATIONS in GFPOP ///////
  ListPiece LP_edges_constraint(Edge const& edge, unsigned int t);
  void LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt);


  ///////  SIMPLE OPERATIONS  ///////
  void addConstant(double myconst);

  void show();

};

#endif // LISTPIECE_H
