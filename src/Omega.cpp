#include "Omega.h"

#include<iostream>
#include <stdlib.h>

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
  LP_ts = NULL;
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Omega::~Omega()
{
  if(LP_ts != NULL)
  {
    for(unsigned int i = 0; i < n + 1; i++)
      {delete[] LP_ts[i]; LP_ts[i] = NULL;}
  }
  delete[] LP_edges;
  LP_edges = NULL;
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
  Interval inter = cost_interval();
  double mini = inter.geta();
  double maxi = inter.getb();
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
    LP_ts[0][i].addFirstPiece(new Piece(Track(), Interval(mini, maxi), Cost()));
    mini = inter.geta();
    maxi = inter.getb();
  }

  ///START STATE CONSTRAINT
  ///START STATE CONSTRAINT
  std::vector<unsigned int> startState = m_graph.getStartState();
  if(startState.size() != 0)
  {
    for(unsigned int i = 0; i < p; i++)
    {
      if(std::find(startState.begin(), startState.end(), i) == startState.end())
        {LP_ts[0][i].setUniquePieceCostToInfinity();}
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
	initialize_LP_ts(n); ///size LP_ts (n+1) x p
	//////////////////////////////

	for(unsigned int t = 0; t < 1; t++) /// loop for all data point (except the first one)
	{
	  LP_edges_operators(t); ///fill_LP_edges. t = newLabel to consider
	  //LP_edges_addPointAndPenalty(myData[t]); ///Add new data point and penalty
	  //LP_t_new_multipleMinimization(t); ///multiple_minimization
	}

	//backtracking();
	show();
}

//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///
//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///

void Omega::LP_edges_operators(unsigned int newLabel)
{
  for(unsigned int i = 0 ; i < q ; i++) /// loop for all edges
  {
    /// i-th edge = m_graph.getEdge(i)
    /// starting state = m_graph.getEdge(i).getState1()
    LP_edges[i].LP_edges_constraint(LP_ts[newLabel][m_graph.getEdge(i).getState1()], m_graph.getEdge(i), newLabel);
  }
}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void Omega::LP_edges_addPointAndPenalty(Point const& pt)
{
  for(unsigned char i = 0; i < q; i++) /// loop for all edges
  {
    /// LP_edges[i] = i-th edge = m_graph.getEdge(i) = we need K, a and penalty
    LP_edges[i].LP_edges_addPointAndPenalty(m_graph.getEdge(i), pt);
  }
}

//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///
//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///

void Omega::LP_t_new_multipleMinimization(unsigned int t)
{
  ///m_graph is rearranged with increasing integer state2
  unsigned int j = 0;
  for (unsigned int i = 0 ; i < p; i++)
  {
    LP_ts[t + 1][i].copy(LP_edges[j]);  //copy pointers in Q_ts[t + 1][i] from Q_edges
    while((j + 1 < q) && (m_graph.getEdge(j + 1).getState2() == i))
    {
      LP_ts[t + 1][i].LP_ts_Minimization(LP_edges[j + 1]); ///
      j = j + 1;
    }
    j = j + 1;
  }
}




//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///
//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///

void Omega::backtracking()
{
  Interval constrainedInterval; ///Interval to fit the constraints
  ///
  /// malsp = Min_Argmin_Label_State_Position
  ///
  double* malsp = LP_ts[n][0].get_min_argmin_label_state_position_ListPiece();
  double* malsp_temp = malsp;

  ///////////////////
  /// FINAL STATE ///
  ///////////////////
  unsigned int CurrentState = 0; ///Current state
  unsigned int CurrentChgpt = n; /// data(1)....data(n). Last data index in each segment
  std::vector<unsigned int> endState = m_graph.getEndState();


  /// IF no endState, all the states are endstates.
  if(endState.size() == 0)
  {
    for (unsigned int j = 1 ; j < p ; j++) ///for all states
    {
      malsp_temp = LP_ts[n][j].get_min_argmin_label_state_position_ListPiece();
      if(malsp_temp[0] < malsp[0]){CurrentState = j; malsp[0] = malsp_temp[0];}
    }
  }
  else
  {
    for (unsigned int j = 0 ; j < endState.size() ; j++) ///for all states
    {
      malsp_temp = LP_ts[n][endState[j]].get_min_argmin_label_state_position_ListPiece();
      if(malsp_temp[0] < malsp[0]){CurrentState = endState[j]; malsp[0] = malsp_temp[0];}
    }
  }

  malsp = LP_ts[n][CurrentState].get_min_argmin_label_state_position_ListPiece();
  globalCost = malsp[0];

  parameters.push_back(malsp[1]);
  changepoints.push_back(CurrentChgpt);
  states.push_back(CurrentState);


  /// BACKTRACK
  ///////////////////////////////
  /// previous to FINAL STATE ///
  ///////////////////////////////

  bool boolForced = false;
  double decay = 0;
  double correction = 1;

  while(malsp[2] > 0) ///while Label > 0
  {
    ///
    ///BACKTRACK
    ///
    boolForced = false;
    decay = m_graph.recursiveState(CurrentState);

    if(decay != 1){correction = std::pow(decay, parameters.back() - malsp[2] + 1);}

    CurrentState = malsp[3];
    CurrentChgpt = malsp[2];

    //TO UPDATE
    //malsp = LP_ts[(int) malsp[2]][(int) malsp[3]].get_min_argmin_label_state_position((int) malsp[4], constrainedInterval, boolForced); ///update boolForced

    //if(malsp[1] > m_bound.getM()){malsp[1] = m_bound.getM(); boolForced = true;}
    //if(malsp[1] < m_bound.getm()){malsp[1] = m_bound.getm(); boolForced = true;}

    parameters.push_back(malsp[1]);
    changepoints.push_back(CurrentChgpt);
    states.push_back(CurrentState);
    forced.push_back(boolForced);
  }


  delete(malsp);
  delete(malsp_temp);
}

///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###


void Omega::show()
{
  for(unsigned char i = 0; i < q; i++)
  {
    std::cout << "s1: " << m_graph.getEdge(i).getState1() + 1;
    std::cout << " s2: " << m_graph.getEdge(i).getState2() + 1 << " ";
    LP_edges[i].show();
  }
}



