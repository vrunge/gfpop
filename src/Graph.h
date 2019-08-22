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

    unsigned int nb_states() const;
    unsigned int nb_edges() const;
    Edge getEdge(unsigned int i) const;
    std::vector<unsigned int> getStartState() const;
    std::vector<unsigned int> getEndState() const;

    bool AreVerticesCompatible() const;
    std::string getType() const;

    Interval buildInterval(double argmin, unsigned int s1, unsigned int s2, bool& out) const;
    double stateDecay(unsigned int s) const;

    void show() const;

    void operator<<(Edge const& newEdge);


  private:
    std::vector<Edge> edges; ///vector edges
    std::vector<unsigned int> startState;
    std::vector<unsigned int> endState;
};

#endif // GRAPH_H
