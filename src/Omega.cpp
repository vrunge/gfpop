#include "Omega.h"
#include "Piece.h"

#include<iostream>
#include <iomanip> ///Set Precision scientific

#include<typeinfo>
#include <stdlib.h>

#include <fstream> ///write in a file
#include<sstream> ///for conversion int to string : function ostringstream

#include <algorithm>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Omega::Omega(Graph graph)
{
  m_graph = graph;
	p = graph.nb_states();
	q = graph.nb_edges();

  /// INITIALIZE ListPiece ///
  LP_edges = new ListPiece[q];
  LP_s_temp = new ListPiece[p];
  LP_ts = NULL;
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Omega::~Omega()
{
  //if(LP_ts != NULL){delete [] LP_ts; LP_ts = NULL;}
  //delete LP_edges;
  //LP_edges = NULL;
  //delete LP_s_temp;
  //LP_s_temp = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > Omega::GetChangepoints() const{return(changepoints);}
std::vector< double > Omega::GetParameters() const{return(parameters);}
std::vector< int > Omega::GetStates() const{return(states);}
std::vector< int > Omega::GetForced() const{return(forced);}
double Omega::GetGlobalCost() const{return(globalCost);}


//####### initialize_LP_ts #######// //####### initialize_LP_ts #######// //####### initialize_LP_ts #######//
//####### initialize_LP_ts #######// //####### initialize_LP_ts #######// //####### initialize_LP_ts #######//

void Omega::initialize_LP_ts(unsigned int n)
{
  double mini = -INFINITY;
  double maxi = INFINITY;
  unsigned int nbR = m_graph.nb_rows();

  LP_ts = new ListPiece*[n+1];
  for(unsigned int i = 0; i < (n+1); i++){LP_ts[i] = new ListPiece[p];}

  ///REVEAL NODE BOUNDARIES IF ANY
  ///REVEAL NODE BOUNDARIES IF ANY
  for(unsigned char i = 0; i < p; i++)
  {
    for(unsigned char j = q; j < nbR; j++)
    {
      if((m_graph.getEdge(j).getConstraint() == "node") && (m_graph.getEdge(j).getState1() == i))
      {
        mini = m_graph.getEdge(j).getMinn();
        maxi = m_graph.getEdge(j).getMaxx();
      }
    }
    LP_ts[0][i].addPiece(new Piece(Track(), Interval(mini, maxi), Cost()));
    mini = -INFINITY;
    maxi = INFINITY;
  }

  ///START STATE CONSTRAINT
  ///START STATE CONSTRAINT
  std::vector<unsigned int> startState = m_graph.getStartState();
  if(startState.size() != 0)
  {
    for(unsigned int i = 0; i < p; i++)
    {
      if(std::find(startState.begin(), startState.end(), i) == startState.end())
        {LP_ts[0][i].addConstant(INFINITY);}
    }
  }
}

//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//

void Omega::gfpop(Data const& data)
{
	Point* myData = data.getVecPt(); ///GET the data = vector of Points = myData
  n = data.getn(); ///data length

	//////////////////////////////
	/// Initialize LP_ts Piece ///
	//////////////////////////////
	initialize_LP_ts(n); ///size LP_ts (n+1) x p

	for(unsigned int t = 0; t < n; t++) /// loop for all data point (except the first one)
	{
	  LP_edges_operators(t); ///fill_LP_edges. t = newLabel to consider
	  LP_edges_addPointAndPenalty(myData[t]); ///Add new data point and penalty
	  //LP_t_new_multipleMinimization(t); ///multiple_minimization
	}

	//backtracking();

	show();
}



//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///

//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///
//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///

void Omega::LP_edges_operators(unsigned int t)
{
  //delete(LP_edges); /// DELETE LP_edges
  unsigned int s1; /// starting state
  Edge edge;
  for(unsigned int i = 0 ; i < q ; i++) /// loop for all edges
  {
    edge = m_graph.getEdge(i);
    s1 = edge.getState1(); /// starting state

    LP_edges[i] = LP_ts[t][s1].LP_edges_constraint(edge, t);
  }
}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void Omega::LP_edges_addPointAndPenalty(Point const& pt)
{
  for(unsigned char i = 0; i < q; i++)
  {
    LP_edges[i].LP_edges_addPointAndPenalty(pt, m_graph.getEdge(i));
  }
}

//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///
//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///

void Omega::LP_t_new_multipleMinimization(unsigned int t)
{
}



void Omega::backtracking()
{

}

///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###


void Omega::show()
{
  for(unsigned char i = 0; i < q; i++){LP_edges[i].show();}
}



