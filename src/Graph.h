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
    std::vector<int> getStartState() const;
    std::vector<int> getEndState() const;

    bool AreVerticesCompatible() const;
    std::string getType() const;

    Interval buildInterval(double argmin, int s1, int s2, bool& out) const;
    double stateDecay(int s) const;

    void show() const;

    void operator<<(Edge const& newEdge);


  private:
    std::vector<Edge> edges; ///vector edges
    std::vector<int> startState;
    std::vector<int> endState;
};

#endif // GRAPH_H
