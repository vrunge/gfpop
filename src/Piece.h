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
    Piece();
    Piece(Track const& info, Interval const& inter = Interval(), Cost const& cost = Cost());
    Piece(const Piece* piece); ///COPY CONSTRUCTOR => copy only the first Piece. piece -> nxt = NULL

    ~Piece();

    Track getTrack() const;
    Interval getInterval() const;
    Cost getCost() const;
    Cost& getRefCost();

    void setIntervalA(double a);
    void setIntervalB(double b);

    Piece* copy();
    double getMin();
    void addCoeff(Cost const& cost, double penalty);

    ///
    ///

    void show();

    Piece *nxt;   /// pointer to next piece

  private:
    Track m_info;
    Interval m_interval;
    Cost m_cost;  /// pointer to the cost associated to the current piece


};

std::ostream &operator>>(std::ostream &flux, Piece* piece);


#endif // PIECE_H
