#ifndef OMEGA_H
#define OMEGA_H

#include"Data.h"
#include"Graph.h"
#include"Edge.h"
#include"Bound.h"
#include"Robust.h"

#include"Piece.h"

#include <math.h>
#include<vector>
#include<list>

#include <stdlib.h>

class Omega
{
  public:
    Omega(Graph graph, Bound bound, Robust robust = Robust());
    ~Omega();

    std::vector< int > GetChangepoints() const;
    std::vector< double > GetMeans() const;
    std::vector< int > GetStates() const;
    std::vector< int > GetForced() const;
    int GetN() const;
    double GetGlobalCost() const;

    void pava(Data const& data);

    void fpop1d_graph_complex(Data const& data);
    void fpop1d_graph_isotonic(Data const& data);
    void fpop1d_graph_std(Data const& data);

    void copyQt();
    void fillQ_edges(int newLabel);
    void multiple_minimization();
    void addPointQ_t(Point pt);

    void backtracking();
    void backtrackingIsotonic(std::vector<Piece*> const& Q_t);

    void save_Q_ts_Q_edges(int t) const;
    void save_Q_s_temp_Q_ts(int t) const;

  private:

    Graph m_graph; ///graph of the constraints
    Bound m_bound; ///min and max of the theta interval. + bool isConstrained = false if all data in [m,M]
    Robust m_robust; ///parameter K and a to define robust loss of type Huber, biweight, L1

    int p;   ///number of states in the graph = number of columns in the matrix Q_ts
    int q; ///number of edges in the graph = number of elements in the object Q_edges

    std::vector< Piece** > Q_ts;  ///cost function Q with respect to position t and state s (size t x p), t = vector size.
    Piece** Q_edges; /// transformed cost by the operators for each edge (size 1 x q)
    Piece** Q_s_temp; /// cost to compare to Q_ts last element of the vector (size 1 x p)

    std::vector< int > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
    std::vector< double > means; ///vector of means build by fpop. size c
    std::vector< int > states; ///vector of states build by fpop. size c
    std::vector< int > forced; ///vector of forced = 0 or 1. 1 = forced value. size c-1
    double globalCost;
};


std::ostream &operator<<(std::ostream &s, const Omega &om);


#endif // OMEGA_H
