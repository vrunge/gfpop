#ifndef OMEGA_H
#define OMEGA_H

#include"Data.h"
#include"Graph.h"
#include"Edge.h"
#include "ListPiece.h"
#include "Piece.h"

#include <math.h>
#include<vector>

#include <stdlib.h>

class Omega
{
  public:
    Omega(Graph graph);
    ~Omega();

    std::vector< std::vector< int > > GetChangepoints() const;
    std::vector< std::vector< double > > GetParameters() const;
    std::vector< std::vector< int > > GetStates() const;
    std::vector< std::vector< int > > GetForced() const;
    std::vector< double > GetGlobalCost() const;

    ///////////////
    void initialize_LP_ts(Point firstData, unsigned int n);
    void gfpop(Data const& data);
    void gfpopTestMode(Data const& data);

    ///////////////
    void LP_edges_operators(unsigned int t);
    void LP_edges_addPointAndPenalty(Point const& pt);
    void LP_t_new_multipleMinimization(unsigned int t);
    void backtracking();
    void show();

  private:
    Graph m_graph; ///graph of the constraints. 9 variables
    unsigned int p;   ///number of states in the graph = number of columns in the matrix Q_ts
    unsigned int q; ///number of edges in the graph = number of elements in the object Q_edges

    unsigned int n; //size of the data
    ListPiece* LP_edges; /// transformed cost by the operators for each edge (size 1 x q)
    ListPiece** LP_ts;  ///cost function Q with respect to position t and state s (size t x p), t = vector size.

    std::vector< std::vector< int > > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
    std::vector< std::vector< double > > parameters; ///vector of means build by fpop. size c
    std::vector< std::vector< int > > states; ///vector of states build by fpop. size c
    std::vector< std::vector< int > > forced; ///vector of forced = 0 or 1. 1 = forced value. size c-1
    std::vector< double > globalCost;
};

std::ostream &operator<<(std::ostream &s, const Omega &om);


#endif // OMEGA_H
