#ifndef OMEGA_H
#define OMEGA_H

#include"Data.h"
#include"Graph.h"
#include"Edge.h"

#include "ListPiece.h"
#include"Piece.h"

#include"Bound.h"
#include"Robust.h"


#include <math.h>
#include<vector>
#include<list>

#include <stdlib.h>

class Omega
{
  public:
    Omega(Graph graph);
    ~Omega();

    std::vector< int > GetChangepoints() const;
    std::vector< double > GetParameters() const;
    std::vector< int > GetStates() const;
    std::vector< int > GetForced() const;
    double GetGlobalCost() const;

    ///////////////
    ///////////////
    ///////////////

    void initialize_LP_ts(unsigned int n);
    void gfpop(Data const& data);

    ///////////////
    ///////////////
    void fill_LP_edges(int newLabel);
    void addPointAndPenalty_LP_t(Point pt);
    void multipleMinimization_LP_edges(int t);

    void show();
    ///////////////
    ///////////////

    void fillQ_edges(int newLabel);
    void multiple_minimization(int t);
    void addPointQ_t(Point pt, int t);

    void backtracking();

  private:
    Graph m_graph; ///graph of the constraints. 9 variables
    unsigned int p;   ///number of states in the graph = number of columns in the matrix Q_ts
    unsigned int q; ///number of edges in the graph = number of elements in the object Q_edges

    unsigned int n; //size of the data
    ListPiece* LP_edges; /// transformed cost by the operators for each edge (size 1 x q)
    ListPiece* LP_s_temp; /// cost to compare to Q_ts last element of the vector (size 1 x p)
    ListPiece** LP_ts;  ///cost function Q with respect to position t and state s (size t x p), t = vector size.

    std::vector< int > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
    std::vector< double > parameters; ///vector of means build by fpop. size c
    std::vector< int > states; ///vector of states build by fpop. size c
    std::vector< int > forced; ///vector of forced = 0 or 1. 1 = forced value. size c-1
    double globalCost;

    ///////////////////////////////////////
    ///////////////TO DELETE///////////////
    ///////////////////////////////////////
    Piece*** Q_ts;  ///cost function Q with respect to position t and state s (size t x p), t = vector size.
    Piece** Q_edges; /// transformed cost by the operators for each edge (size 1 x q)
    Piece** Q_s_temp; /// cost to compare to Q_ts last element of the vector (size 1 x p)

    Bound m_bound; ///min and max of the theta interval. + bool isConstrained = false if all data in [m,M]
    Robust m_robust; ///parameter K and a to define robust loss of type Huber, biweight, L1
};


std::ostream &operator<<(std::ostream &s, const Omega &om);


#endif // OMEGA_H
