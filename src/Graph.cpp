#include "Graph.h"

#include "termcolor.h"

#include<iostream>
#include<string>
#include<vector>
#include <math.h>

#include <algorithm>

Graph::Graph()
{
  startState = -1; ///undefined start state
  endState = -1; ///undefined end state
}

void Graph::newEdge(Edge const& edge){edges.push_back(edge);}

// ### nb_states ### ////// ### nb_states ### ////// ### nb_states ### ////// ### nb_states ### ///
// ### nb_states ### ////// ### nb_states ### ////// ### nb_states ### ////// ### nb_states ### ///

int Graph::nb_states() const
{
  std::vector<int> temp;
  for (int i = 0 ; i < edges.size() ; i++)
  {
    temp.push_back(edges[i].getState1());
    temp.push_back(edges[i].getState2());
  }
  sort(temp.begin(), temp.end());

  int res = 1;
  for(int j = 0; j < temp.size() - 1 ; j++)
  {
    if(temp[j] != temp[j + 1]){res = res + 1;}
  }
  return(res);
}


// ### get ### ////// ### get ### ////// ### get ### ////// ### get ### ///
// ### get ### ////// ### get ### ////// ### get ### ////// ### get ### ///

int Graph::nb_edges() const {return(edges.size());}
Edge Graph::getEdge(int i) const {return(edges[i]);}
int Graph::getStartState() const {return(startState);}
int Graph::getEndState() const {return(endState);}


// ### isCyclic ### ////// ### isCyclic ### ////// ### isCyclic ### ////// ### isCyclic ### ///
// ### isCyclic ### ////// ### isCyclic ### ////// ### isCyclic ### ////// ### isCyclic ### ///

bool Graph::isCyclic() const
{
  int nbEdges = nb_edges();
  if(nbEdges == 1){return(false);}
  if(nbEdges!= nb_states()){return(false);}
  std::vector<int> state1(nbEdges, 0);
  std::vector<int> state2(nbEdges, 0);

  for (int i = 0 ; i < nbEdges ; i++)
  {
    if(state1[edges[i].getState1()] == 1){return(false);}
    if(state2[edges[i].getState2()] == 1){return(false);}
    state1[edges[i].getState1()] = 1;
    state2[edges[i].getState2()] = 1;
  }
  return(true);
}


bool Graph::AreVerticesCompatible() const //label of the vertices from 0 to S
{
  int maxLabel = 0;
  int nbEdges = nb_edges();
  for (int i = 0 ; i < nbEdges ; i++)
  {
    if(edges[i].getState1() > maxLabel){maxLabel = edges[i].getState1();}
    if(edges[i].getState2() > maxLabel){maxLabel = edges[i].getState2();}
  }

  bool isAvertex[maxLabel + 1] = {false};
  for (int i = 0 ; i < nbEdges ; i++)
  {
    isAvertex[edges[i].getState1()] = true;
    isAvertex[edges[i].getState2()] = true;
  }

  int nb = 0;
  int j = 0;
  while((j <= maxLabel) && (isAvertex[j] == true)){nb = nb + 1; j = j + 1; }
  return(nb == nb_states());
}


// ### getType ### /// /// ### getType ### /// /// ### getType ### /// /// ### getType ### ///
// ### getType ### /// /// ### getType ### /// /// ### getType ### /// /// ### getType ### ///


std::string Graph::getType() const
{
  std::string response = "complex";
  if(edges.size() == 1)
  {
    if(edges[0].getConstraint() == "up"){response = "isotonic";}
    if(edges[0].getConstraint() == "std" && nb_states() == 1){response = "std";}
  }

  if(isCyclic() == true){response = "cyclic";}

  return(response);

}


// ### buildInterval ### ////// ### buildInterval ### ////// ### buildInterval ### ////// ### buildInterval ### ///
// ### buildInterval ### ////// ### buildInterval ### ////// ### buildInterval ### ////// ### buildInterval ### ///


Interval Graph::buildInterval(double argmin, int s1, int s2, bool& out) const
{
  Interval response = Interval(-INFINITY, INFINITY);

  /// FIND edge
  Edge myedge;
  for (int i = 0 ; i < edges.size() ; i++)
  {
    if((edges[i].getState1() == s1) && (edges[i].getState2() == s2)){myedge = edges[i];}
  }

  if(myedge.getConstraint() == "up")
  {
    response.setb(argmin - myedge.getParameter());
  }

  if(myedge.getConstraint()  == "down")
  {
    response.seta(argmin + myedge.getParameter());
  }

  if(myedge.getConstraint()  == "absInf")
  {
    response.seta(argmin - myedge.getParameter());
    response.setb(argmin + myedge.getParameter());
  }

  if(myedge.getConstraint()  == "absSup") ///DANGER : exclusion of this interval
  {
    response.seta(argmin - myedge.getParameter());
    response.setb(argmin + myedge.getParameter());
    out = true;
  }

  return(response);
}


// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///
// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///


void Graph::show() const
{
  //std::cout<< termcolor::on_red << "GRAPH" << termcolor::reset << std::endl;
  for (int i = 0 ; i < edges.size() ; i++)
  {
    //edges[i].show();
  }
  //std::cout<< "Start state : " << startState << std::endl;
  //std::cout<< "End state : " << endState << std::endl;
  //std::cout<< "nb states : " << nb_states() << std::endl;
  //std::cout<< "nb edges : " << nb_edges() << std::endl;
}

// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///
// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///

void Graph::operator<<(Edge const& newEdge)
{
  if(newEdge.getConstraint() == "start"){startState = newEdge.getState1();}
  if(newEdge.getConstraint() == "end"){endState = newEdge.getState1();}
  if((newEdge.getConstraint() != "start") && (newEdge.getConstraint() != "end")){edges.push_back(newEdge);}
}
