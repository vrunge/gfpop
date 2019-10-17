/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "Graph.h"

#include<iostream>
#include<string>
#include<vector>
#include <math.h>

#include <algorithm>

Graph::Graph()
{
}

void Graph::newEdge(Edge const& edge){edges.push_back(edge);}

// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### ///
// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### ///

unsigned int Graph::nb_states() const
{
  std::vector<unsigned int> temp;
  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    temp.push_back(it->getState1());
    temp.push_back(it->getState2());
  }
  sort(temp.begin(), temp.end());

  unsigned int res = 1;
  for(unsigned int j = 0; j < temp.size() - 1 ; j++)
  {
    if(temp[j] != temp[j + 1]){res = res + 1;}
  }
  return(res);
}


// ### nb_edges ### /// /// ### nb_edges ### /// /// ### nb_edges ### /// /// ### nb_edges ### ///
// ### nb_edges ### /// /// ### nb_edges ### /// /// ### nb_edges ### /// /// ### nb_edges ### ///

unsigned int Graph::nb_edges() const
{
  unsigned int res = 0;
  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    if(edges[i].getConstraint() != "node"){res = res + 1;};
  }
  return(res);
}


// ### nb_rows ### /// /// ### nb_rows ### /// /// ### nb_rows ### ////// ### nb_rows ### ///
// ### nb_rows ### /// /// ### nb_rows ### /// /// ### nb_rows ### ////// ### nb_rows ### ///

unsigned int Graph::nb_rows() const {return(edges.size());}


// ### get ### /// /// ### get ### /// /// ### get ### ////// ### get ### ///
// ### get ### /// /// ### get ### /// /// ### get ### ////// ### get ### ///

Edge Graph::getEdge(unsigned int i) const {return(edges[i]);}
std::vector<unsigned int> Graph::getStartState() const {return(startState);}
std::vector<unsigned int> Graph::getEndState() const {return(endState);}


// ### AreVerticesCompatible ### /// /// ### AreVerticesCompatible ### /// /// ### AreVerticesCompatible ### ///
// ### AreVerticesCompatible ### /// /// ### AreVerticesCompatible ### /// /// ### AreVerticesCompatible ### ///


bool Graph::AreVerticesCompatible() const //label of the vertices from 0 to S
{
  unsigned int maxLabel = 0;
  unsigned int nbEdges = nb_edges();
  for (unsigned int i = 0 ; i < nbEdges ; i++)
  {
    if(edges[i].getState1() > maxLabel){maxLabel = edges[i].getState1();}
    if(edges[i].getState2() > maxLabel){maxLabel = edges[i].getState2();}
  }

  bool isAvertex[maxLabel + 1];
  for(unsigned int i = 0; i < maxLabel + 1; i++){isAvertex[i] = false;}

  for (unsigned int i = 0 ; i < nbEdges ; i++)
  {
    isAvertex[edges[i].getState1()] = true;
    isAvertex[edges[i].getState2()] = true;
  }

  unsigned int nb = 0;
  unsigned int j = 0;
  while((j <= maxLabel) && (isAvertex[j] == true)){nb = nb + 1; j = j + 1;}
  return(nb == nb_states());
}


// ### getType ### /// /// ### getType ### /// /// ### getType ### /// /// ### getType ### ///
// ### getType ### /// /// ### getType ### /// /// ### getType ### /// /// ### getType ### ///


std::string Graph::getType() const
{
  std::string response = "complex";
  if(edges.size() == 2)
  {
    if(edges[0].getConstraint() == "null" && edges[1].getConstraint() == "up" && nb_states() == 1 && edges[0].getParameter() == 1){response = "isotonic";}
    if(edges[0].getConstraint() == "null" && edges[1].getConstraint() == "up" && nb_states() == 1 && edges[0].getParameter() == 1  && edges[1].getBeta() == 0){response = "pava";}
    if(edges[0].getConstraint() == "null" && edges[1].getConstraint() == "std" && nb_states() == 1 && edges[0].getParameter() == 1){response = "std";}
   }
  return(response);
}


// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///
// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///


Interval Graph::buildInterval(double argmin, unsigned int s1, unsigned int s2, bool& out) const
{
  Interval response = Interval(-INFINITY, INFINITY);

  /// FIND edge. If there are 2 edges (s1,s2) we get the second one (which is not of "null" or "decay" type). (cf ordering in gfpop R function)
  Edge myedge;
  for (unsigned int i = 0 ; i < edges.size() ; i++)
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



// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///
// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///


double Graph::stateDecay(unsigned int s) const
{
  double response =  1;
  for (int i = edges.size()-1 ; i > -1; i--)
  {
    if((edges[i].getState1() == s) && (edges[i].getState2() == s)){response = edges[i].getParameter();}
  }

  return(response);
}

// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///
// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///


void Graph::show() const
{
  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    edges[i].show();
  }
  for (unsigned int i = 0 ; i < startState.size() ; i++)
  {
    std::cout<< startState[i] << " ";
  }
  for (unsigned int i = 0 ; i < endState.size() ; i++)
  {
    std::cout<< endState[i] << " ";
  }
  std::cout << std::endl;
  std::cout<< "nb states : " << nb_states() << std::endl;
  std::cout<< "nb edges : " << nb_edges() << std::endl;
}

// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///
// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///

void Graph::operator<<(Edge const& newEdge)
{
  if(newEdge.getConstraint() == "start"){startState.push_back(newEdge.getState1());}
  if(newEdge.getConstraint() == "end"){endState.push_back(newEdge.getState1());}
  if((newEdge.getConstraint() != "start") && (newEdge.getConstraint() != "end")){edges.push_back(newEdge);}
}


