/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "Graph.h"

#include<iostream>

Graph::Graph(){}

void Graph::newEdge(Edge const& edge){edges.push_back(edge);}

// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### ///
// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### /// /// ### nb_states ### ///

unsigned int Graph::nb_states() const
{
  std::vector<unsigned int> temp;
  for (std::vector<Edge>::const_iterator it=edges.begin(); it < edges.end(); it++)
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


// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///
// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### /// /// ### buildInterval ### ///

Interval Graph::buildInterval(double argmin, unsigned int s1, unsigned int s2, bool& out) const
{
  Interval response = cost_interval();
  unsigned int nb = 0;
  unsigned int edgeIndex;
  Interval inter = cost_interval();

  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    if((edges[i].getState1() == s1) && (edges[i].getState2() == s2))
    {
      cost_interShift(argmin, - edges[i].getParameter());
      if(edges[i].getConstraint() == "up"){response.setb(cost_interShift(argmin, -edges[i].getParameter())); nb = nb + 1; edgeIndex = i;}
      if(edges[i].getConstraint()  == "down"){response.seta(cost_interShift(argmin, edges[i].getParameter())); nb = nb + 1;  edgeIndex = i;}
      if(edges[i].getConstraint()  == "node"){inter = Interval(edges[i].getMinn(), edges[i].getMaxx());}
    }
  }

  if(nb == 2) /// abs (= up + down edges) case
  {
    out = true;
    response.seta(cost_interShift(argmin, - edges[edgeIndex].getParameter()));
    response.setb(cost_interShift(argmin, edges[edgeIndex].getParameter()));
  }

  response.seta(std::max(inter.geta(), response.geta()));
  response.setb(std::min(inter.getb(), response.getb()));

  return(response);
}


// ### recursiveState ### /// /// ### recursiveState ### /// /// ### recursiveState ### /// /// ### recursiveState ### ///
// ### recursiveState ### /// /// ### recursiveState ### /// /// ### recursiveState ### /// /// ### recursiveState ### ///

double Graph::recursiveState(unsigned int s) const
{
  double response =  1;
  for(unsigned int i = 0; i < edges.size(); i++)
  {
    if((edges[i].getState1() == s) && (edges[i].getState2() == s) && (edges[i].getConstraint()  == "null")){response = edges[i].getParameter();}
  }
  return(response);
}


// ### findBeta ### /// /// ### findBeta ### /// /// ### findBeta ### /// /// ### findBeta ### ///
// ### findBeta ### /// /// ### findBeta ### /// /// ### findBeta ### /// /// ### findBeta ### ///

double Graph::findBeta(unsigned int state1, unsigned int state2)
{
  double response = 0;
  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    if((edges[i].getState1() == state1) && (edges[i].getState2() == state2) && (edges[i].getConstraint() != "node"))
    {
      response = edges[i].getBeta();
    }
  }
  return(response);
}



// ### nodeConstraints ### /// /// ### nodeConstraints ### /// /// ### nodeConstraints ### /// /// ### nodeConstraints ### ///
// ### nodeConstraints ### /// /// ### nodeConstraints ### /// /// ### nodeConstraints ### /// /// ### nodeConstraints ### ///

Interval* Graph::nodeConstraints()
{
  Interval* inter = new Interval[nb_states()];
  for (unsigned int i = 0 ; i < nb_states(); i++)
  {
    inter[i] = cost_interval();
  }
  for (unsigned int i = 0 ; i < edges.size(); i++)
  {
    if(edges[i].getConstraint() == "node"){inter[edges[i].getState1()] = Interval(edges[i].getMinn(), edges[i].getMaxx());}
  }
  return(inter);
}


// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///
// ### show ### /// /// ### show ### /// /// ### show ### /// /// ### show ### ///

void Graph::show() const
{
  for (unsigned int i = 0 ; i < edges.size() ; i++)
  {
    edges[i].show();
  }
  //std::cout << "Start state (+1): ";
  for (unsigned int i = 0 ; i < startState.size() ; i++)
  {
    //std::cout << startState[i] << " ";
  }
  //std::cout << std::endl;
  //std::cout<< "End state (+1): ";
  for (unsigned int i = 0 ; i < endState.size() ; i++)
  {
    //std::cout << endState[i] << " ";
  }
  //std::cout << std::endl;
  //std::cout<< "nb states : " << nb_states() << std::endl;
  //std::cout<< "nb edges : " << nb_edges() << std::endl;
}

// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///
// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### /// /// ### OPERATOR ### ///

void Graph::operator<<(Edge const& newEdge)
{
  if(newEdge.getConstraint() == "start"){startState.push_back(newEdge.getState1());}
  if(newEdge.getConstraint() == "end"){endState.push_back(newEdge.getState1());}
  if((newEdge.getConstraint() != "start") && (newEdge.getConstraint() != "end")){edges.push_back(newEdge);}
}


