#ifndef PIECE_H
#define PIECE_H

#include "Edge.h"
#include "Data.h"

#include "Track.h"
#include "Interval.h"
#include "Cost.h"

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
    double getMin();

    void addCostAndPenalty(Cost const& cost, double penalty);

    ///
    ///
    Interval intervalMinLessUp(double bound, double currentValue, bool constPiece);
    Interval intervalMinLessDw(double bound, double currentValue, bool constPiece);
    Piece* pastePieceUp(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack);
    Piece* pastePieceDw(const Piece* NXTPiece, Interval const& decrInter, Track const& newTrack);
    ///
    ///

    void show();


};

std::ostream &operator>>(std::ostream &flux, Piece* piece);


#endif // PIECE_H
