#ifndef PIECE_H
#define PIECE_H

#include "Edge.h"
#include "Bound.h"
#include "Robust.h"

#include "Track.h"
#include "Interval.h"
#include "CostGauss.h"

#include<vector>
#include<string>

#include <fstream> ///write in a file

class Piece
{
  public:
    Piece();
    Piece(Track const& info, Interval const& inter = Interval(), CostGauss const& cost = CostGauss());
    Piece(const Piece* piece); ///COPY CONSTRUCTOR => copy only the first Piece. piece -> nxt = NULL

    ~Piece();

    Track getTrack() const;
    Interval getInterval() const;
    CostGauss getCost() const;
    CostGauss& getRefCost();

    int length();
    //###################


    void testEmptyInterval();

    Piece* copy();
    Piece* copy(int& length);
    Piece* copyIsotonic(double newLeftBound);
    Piece* reverse();
    Piece* reverse(int length);
    void opposition();

   //###################

    double* Piece_min_argmin();


    double newBound(double newleftBound);

    //###################

    void addConstant(double cst);
    void addPoint(Point const& pt, Robust const& robust);

    //###################

    Piece* edge_constraint(Edge const& edge, int newLabel, Bound const& bound);

    //###################

    Piece* shift_right(double parameter, double maxi);
    Piece* shift_left(double parameter, double mini);

      Piece* operator_std_min_argmin(int newLabel, int& labelmin, double& argmini, Bound const& bound);
      Piece* operator_stdConstr_min_argmin(int newLabel, int& labelmin, double& argmini, Bound const& bound);

      Piece* operator_std(int newLabel, int parentStateLabel, Bound const& bound);
      Piece* operator_stdConstr(int newLabel, int parentStateLabel, Bound const& bound);
      Piece* operator_up(int newLabel, int edgeLabel);
      Piece* operator_down(int newLabel, int edgeLabel);
      Piece* pastePiece(const Piece* Q, Interval const& decrInter, Track const& newTrack);

    //###################

    Piece* min_function(Piece* Q2, double M);
    Piece* max_function(Piece* Q2, double M);
    Piece* pieceGenerator(Piece* Q1, Piece* Q2, int Bound_Q2_Minus_Q1, double M);

    double* get_min_argmin_label_state_position_final();
    double* get_min_argmin_label_state_position(int i, Interval const& constrainedInter, bool out, bool& forced, bool isBoundConstrained);
    double* get_min_argmin_label(double rightBound, bool& forced, bool isBoundConstrained);
    //###################

    void save(std::ostream &flux);
    void show();
    void showOne();

  private:

    Track m_info;

    Interval m_interval;
    CostGauss m_cost;  /// pointer to the cost associated to the current piece

    Piece *nxt;   /// pointer to next piece

};

std::ostream &operator>>(std::ostream &flux, Piece* piece);


#endif // PIECE_H
