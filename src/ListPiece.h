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
  void setUniquePieceCostToInfinity();
  void reset();
  void reverse();
  void move();
  void initializeCurrentPiece();
  void addCurrentPiecePlus1(Piece* newPiece);
  void addNewLastPiece(Piece* newPiece);
  void copy(ListPiece  const& LP_edge);

  void shift(double parameter);
  void expDecay(double gamma);

  ///////  3 OPERATIONS in GFPOP ///////
  void LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel);
  void LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt);
  void LP_ts_Minimization(ListPiece const& LP_edge);

  void operator_up(ListPiece const& LP_edge, unsigned int newLabel, unsigned int parentState, bool upDirection);

  void show();

};

#endif // LISTPIECE_H
