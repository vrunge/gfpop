//  GPL-3 License
// Copyright (c) 2019 Vincent Runge

#ifndef GRAPH_H
#define GRAPH_H

#include"Edge.h"
#include"Interval.h"
#include"ExternFunctions.h"

#include<vector>
#include<string>
#include <math.h>
#include <algorithm>


class Graph
{
  public:
    Graph();
    void newEdge(Edge const& edge);

    unsigned int nb_states() const;
    unsigned int nb_edges() const;
    unsigned int nb_rows() const;

    Edge getEdge(unsigned int i) const;
    std::vector<unsigned int> getStartState() const;
    std::vector<unsigned int> getEndState() const;

    Interval buildInterval(double argmin, unsigned int s1, unsigned int s2, bool& out) const;
    double recursiveState(unsigned int s) const;
    double findBeta(unsigned int state1, unsigned int state2);
    Interval* nodeConstraints();

    void show() const;

    void operator<<(Edge const& newEdge);


  private:
    std::vector<Edge> edges; ///vector edges
    std::vector<unsigned int> startState;
    std::vector<unsigned int> endState;
};

#endif // GRAPH_H
