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
  void setNewBounds(Interval newBounds);

  void reset();
  void copy(ListPiece  const& LP_edge);

  void reverseAndCount(unsigned int& length);
  void reverseAndSetTrackPosition(unsigned int length);

  void move();
  void initializeCurrentPiece();

  void addCurrentPiecePlus1NotMove(Piece* newPiece);
  void addFirstPiece(Piece* newPiece);

  void shift(double parameter);
  void expDecay(double gamma);

  ///////  3 OPERATIONS in GFPOP ///////
  void LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel);
  void LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt);
  void LP_ts_Minimization(ListPiece& LP_edge);

  ///////  operators up and down ///////
  void operatorUp(ListPiece const& LP_edge, unsigned int newLabel, unsigned int parentState);
  void operatorDw(ListPiece const& LP_edge, unsigned int newLabel, unsigned int parentState);

  ///////  get info ///////
  double* get_min_argmin_label_state_position_ListPiece();
  double* get_min_argmin_label_state_position_onePiece(unsigned int position, Interval constrainedInterval, bool& forced);

  void show();

};

#endif // LISTPIECE_H
