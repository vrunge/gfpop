#ifndef GRAPH_H
#define GRAPH_H

#include"Edge.h"
#include"Interval.h"

#include<vector>
#include<string>

class Graph
{
  public:
    Graph();
    void newEdge(Edge const& edge);

    int nb_states() const;
    int nb_edges() const;
    Edge getEdge(int i) const;
    int getStartState() const;
    int getEndState() const;

    bool isCyclic() const;
    std::string getType() const;

    Interval buildInterval(double argmin, int s1, int s2, bool& out) const;

    void show() const;

    void operator<<(Edge const& newEdge);


  private:
    std::vector<Edge> edges; ///vector edges
    int startState;
    int endState;
};

#endif // GRAPH_H
