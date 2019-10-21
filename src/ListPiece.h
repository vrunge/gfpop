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
  Piece* lastPiece;

public:
  ListPiece();
  ~ListPiece();

  ///////  Simple list operations  ///////
  void reset();
  void reverse();
  void move();
  void initializeCurrentPiece();
  void addCurrentPiecePlus1(Piece* newPiece);
  void addNewLastPiece(Piece* newPiece);
  void copy(ListPiece  const& LP_edge);

  void shift_right(double parameter);
  void shift_left(double parameter);

  ///////  3 OPERATIONS in GFPOP ///////
  void LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel);
  void LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt);
  void LP_ts_Minimization(ListPiece const& LP_edge);


  ///////  Simple Piece operations  ///////
  void setUniquePieceCostToInfinity();

  void show();

};

#endif // LISTPIECE_H
