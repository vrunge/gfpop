#ifndef PIECE_H
#define PIECE_H

#include "Edge.h"
#include "Data.h"

#include "Track.h"
#include "Interval.h"
#include "Cost.h"
#include <math.h>

#include<vector>
#include<string>

#include <fstream> ///write in a file

class Piece
{
  public:
    Track m_info;
    Interval m_interval;
    Cost m_cost;
    Piece *nxt;   /// pointer to next piece

    Piece();
    Piece(Track const& info, Interval const& inter = Interval(), Cost const& cost = Cost());
    Piece(const Piece* piece); ///COPY CONSTRUCTOR => copy only the first Piece. piece -> nxt = NULL
    ~Piece();
    Piece* copy();

    void addCostAndPenalty(Cost const& cost, double penalty);

    ///
    ///
    Interval intervalMinLessUp(double bound, double currentValue, bool constPiece);
    Interval intervalMinLessDw(double bound, double currentValue, bool constPiece);
    Piece* pastePieceUp(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack);
    Piece* pastePieceDw(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack);

    Piece* pieceGenerator(Piece* Q1, Piece* Q2, int Bound_Q2_Minus_Q1, double M);
    Piece* piece0(Piece* Q1, Piece* Q2, Interval interToPaste, int& Q2_Minus_Q1);
    Piece* piece1(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1);
    Piece* piece2(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1);

    void get_min_argmin_label_state_position(double* response);
    ///
    ///

    void show();

};

std::ostream &operator>>(std::ostream &flux, Piece* piece);


#endif // PIECE_H
