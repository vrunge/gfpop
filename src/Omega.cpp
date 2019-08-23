#include "Omega.h"
#include "CostGauss.h"
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
  //////////TO DELETE//////////
  m_bound = Bound(0, 0, false);
  m_robust = Robust(1000,1000);

	p = graph.nb_states();
	q = graph.nb_edges();
  double mini = -INFINITY;
  double maxi = INFINITY;
  unsigned int nbR = graph.nb_rows();

  Q_ts = NULL;

  /// INITIALIZE Q_edges ///
	Q_edges = new Piece*[q];
	for(unsigned char i = 0; i < q; i++){Q_edges[i] = new Piece(Track(), Interval(-INFINITY, INFINITY), CostGauss());}

	/// INITIALIZE Q_s_temp ///
	Q_s_temp = new Piece*[p];
  for(unsigned char i = 0; i < p; i++)
  {
    ///REVEAL NODE BOUNDARIES IF ANY
    for(unsigned char j = q; j < nbR; j++)
    {
      if((graph.getEdge(j).getConstraint() == "node") && (graph.getEdge(j).getState1() == i))
        {
          mini = graph.getEdge(j).getMinn();
          maxi = graph.getEdge(j).getMaxx();
        }
    }

    Q_s_temp[i] = new Piece(Track(), Interval(mini, maxi), CostGauss());
    mini = INFINITY;
    maxi = -INFINITY;
  }
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Omega::~Omega()
{
  if(Q_ts != NULL){for(unsigned int i = 0; i < (n + 1); i++){delete [] Q_ts[i]; Q_ts[i] = NULL;}}
  delete [] Q_edges;
  Q_edges = NULL;
  delete [] Q_s_temp;
  Q_s_temp = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > Omega::GetChangepoints() const{return(changepoints);}
std::vector< double > Omega::GetParameters() const{return(parameters);}
std::vector< int > Omega::GetStates() const{return(states);}
std::vector< int > Omega::GetForced() const{return(forced);}

int Omega::GetN() const{return(n);}
double Omega::GetGlobalCost() const{return(globalCost);}


//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//
//####### gfpop #######// //####### gfpop #######// //####### gfpop #######//

void Omega::gfpop(Data const& data)
{
	Point* myData = data.getVecPt(); ///GET the data/// get the vector of Points = myData
  n = data.getn();

	///
	/// Initialize Q_ts Piece***
	///
	Q_ts = new Piece**[n + 1];
	for(unsigned int i = 0 ; i < (n + 1) ; i++){Q_ts[i] = new Piece*[p];}

	/// Initialize first functional cost. Interval / add first point / constraint by starting vertices
	for (unsigned char i = 0; i < p ; i++){Q_ts[1][i] = Q_s_temp[0] -> copy();}
	addPointQ_t(myData[0], 0);
  std::vector<unsigned int> startState = m_graph.getStartState();
  if(startState.size() != 0){for(unsigned int i = 0; i < p; i++){if(std::find(startState.begin(), startState.end(), i) == startState.end()){Q_ts[1][i] -> addConstant(INFINITY);}}}

  for(unsigned int t = 1; t < n; t++) /// loop for all data point (except the first one)
  {
    fillQ_edges(t); ///fillQ_edges. t = newLabel to consider
    multiple_minimization(t); ///multiple_minimization
    addPointQ_t(myData[t], t); ///Add new data point
  }
  backtracking();
}



//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///

//##### fillQ_edges #####//////##### fillQ_edges #####//////##### fillQ_edges #####///
//##### fillQ_edges #####//////##### fillQ_edges #####//////##### fillQ_edges #####///

void Omega::fillQ_edges(int newLabel)
{
	int s1; /// starting state
	for (unsigned int i = 0 ; i < q ; i++) /// loop for all edges
	{
    delete(Q_edges[i]); /// DELETE Q_edges[i]
		Edge edge = m_graph.getEdge(i);
    s1 = edge.getState1(); /// starting state
    Q_edges[i] = Q_ts[newLabel][s1] -> edge_constraint(edge, newLabel, m_bound);
	}
}

//##### multiple_minimization #####//////##### multiple_minimization #####//////##### multiple_minimization #####///
//##### multiple_minimization #####//////##### multiple_minimization #####//////##### multiple_minimization #####///

void Omega::multiple_minimization(int t)
{
  unsigned int j = 0;
  /// Q_s_temp vs Q_ts minimization
  for (unsigned int i = 0 ; i < p; i++)
  {
    //copy pointers in Q_ts[t + 1][i] from Q_edges
    Q_ts[t + 1][i] = Q_edges[j] -> copy();
    while((j + 1 < q) && (m_graph.getEdge(j + 1).getState2() == i))
    {
      Q_ts[t + 1][i] = Q_ts[t + 1][i] -> min_function(Q_edges[j + 1], m_bound.getM()); ///
      j = j + 1;
    }
    j = j + 1;
  }
}

//##### addPointQ_t #####//////##### addPointQ_t #####//////##### addPointQ_t #####///
//##### addPointQ_t #####//////##### addPointQ_t #####//////##### addPointQ_t #####///

void Omega::addPointQ_t(Point pt, int t)
{
	for(unsigned char i = 0; i < p; i++){Q_ts[t + 1][i] -> addPoint(pt, m_robust);}
}


//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///
//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///

void Omega::backtracking()
{
  Interval constrainedInterval; ///Interval to fit the constraints
  ///
  /// malsp = Min_Argmin_Label_State_Position_Final
  ///
  std::vector<double> malsp = Q_ts[n][0] -> get_min_argmin_label_state_position_final();
  std::vector<double> malsp_temp;

  ///
  ///FINAL STATE
  ///
  int CurrentState = 0; ///Current state
  int CurrentChgpt = n; /// data(1)....data(n). Last data index in each segment
  std::vector<unsigned int> endState = m_graph.getEndState();

  if(endState.size() == 0)
  {
    for (unsigned int j = 1 ; j < p ; j++) ///for all states
    {
      malsp_temp = Q_ts[n][j] -> get_min_argmin_label_state_position_final();
      if(malsp_temp[0] < malsp[0]){CurrentState = j; malsp[0] = malsp_temp[0];}
    }
  }
  else
  {
    for (unsigned int j = 0 ; j < endState.size() ; j++) ///for all states
    {
      malsp_temp = Q_ts[n][endState[j]] -> get_min_argmin_label_state_position_final();
      if(malsp_temp[0] < malsp[0]){CurrentState = endState[j]; malsp[0] = malsp_temp[0];}
    }
  }

  malsp = Q_ts[n][CurrentState] -> get_min_argmin_label_state_position_final();
  globalCost = malsp[0];

  parameters.push_back(malsp[1]);
  changepoints.push_back(CurrentChgpt);
  states.push_back(CurrentState);


  /// BACKTRACK
  ///
  ///BEFORE FINAL STATE
  ///

  bool out = false;
  bool boolForced = false;
  double decay = 1;
  double correction = 1;


  while(malsp[2] > 0) ///while Label > 0
  {
    ///
    ///BACKTRACK
    ///
    out = false;
    boolForced = false;
    decay = m_graph.stateDecay(CurrentState);
    if(decay != 1){correction = std::pow(decay, parameters.back() - malsp[2] + 1);}

    constrainedInterval = m_graph.buildInterval(malsp[1]*correction, malsp[3], CurrentState, out); ///update out
    CurrentState = malsp[3];
    CurrentChgpt = malsp[2];

    malsp = Q_ts[(int) malsp[2]][(int) malsp[3]] -> get_min_argmin_label_state_position((int) malsp[4], constrainedInterval, out, boolForced, m_bound.getIsConstrained()); ///update boolForced

    if(malsp[1] > m_bound.getM()){malsp[1] = m_bound.getM(); boolForced = true;}
    if(malsp[1] < m_bound.getm()){malsp[1] = m_bound.getm(); boolForced = true;}

    parameters.push_back(malsp[1]);
    changepoints.push_back(CurrentChgpt);
    states.push_back(CurrentState);
    forced.push_back(boolForced);
  }
  //std::cout<<"DONE"<<std::endl;
}



///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###

std::ostream &operator<<(std::ostream &s, const Omega &om)
{
  std::vector< int > chpt = om.GetChangepoints();
  std::vector< double > parameters = om.GetParameters();
  std::vector< int > states = om.GetStates();
  std::vector< int > forced = om.GetForced();
  s << " n : " << om.GetN()-1 << std::endl;
  for (int i = chpt.size()-1; i > -1; i--){s << " ** " << chpt[i];}
  s << std::endl;
  for (int i = parameters.size()-1; i > -1; i--){s << " ** " << parameters[i];}
  s << std::endl;
  for (int i = states.size()-1; i > -1; i--){s << " ** " << states[i];}
  s << std::endl;
  for (int i = forced.size()-1; i > -1; i--){s << " ** " << forced[i];}
  s << std::endl;
  s << "globalCost: "<< om.GetGlobalCost() << std::endl;
  return s;
}


