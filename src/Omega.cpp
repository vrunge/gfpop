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
    for(unsigned int i = 0; i < (n + 1); i++){delete[] (LP_ts[i]);}
    delete[] LP_ts;
    LP_ts = NULL;
  }
  delete[] LP_edges;
  LP_edges = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//
std::vector< std::vector< int > > Omega::GetChangepoints() const{return(changepoints);}
std::vector< std::vector< double > > Omega::GetParameters() const{return(parameters);}
std::vector< std::vector< int > > Omega::GetStates() const{return(states);}
std::vector< std::vector< int > > Omega::GetForced() const{return(forced);}
std::vector< double > Omega::GetGlobalCost() const{return(globalCost);}

//####### initialize_LP_ts #######// //####### initialize_LP_ts #######// //####### initialize_LP_ts #######//
//####### initialize_LP_ts #######// //####### initialize_LP_ts #######// //####### initialize_LP_ts #######//
// initialize LP_ts for all t and s
// t = 0 for all s : LP_ts[0][s] = addFirstPiece(new Piece(Track(), Interval(mini, maxi), 0 or +INFINITY));
// t > 1 for all s : LP_ts[t][s] = addFirstPiece(new Piece(Track(), Interval(mini, maxi), +INFINITY));

void Omega::initialize_LP_ts(Point firstData, unsigned int n)
{
  Interval inter = cost_interval(); ///get the cost-dependent interval
  double mini = inter.geta();
  double maxi = inter.getb();
  unsigned int nbR = m_graph.nb_rows();

  LP_ts = new ListPiece*[n + 1];
  for(unsigned int i = 0; i < (n + 1); i++){LP_ts[i] = new ListPiece[p]; for(unsigned int j = 0; j < p; j++){LP_ts[i][j] = ListPiece();}}

  ///REVEAL NODE BOUNDARIES IF ANY
  ///REVEAL NODE BOUNDARIES IF ANY
  for(unsigned char j = 0; j < p; j++)
  {
    for(unsigned char k = q; k < nbR; k++) ///after the q edges
    {
      if((m_graph.getEdge(k).getConstraint() == "node") && (m_graph.getEdge(k).getState1() == j))
      {
        mini = m_graph.getEdge(k).getMinn();
        maxi = m_graph.getEdge(k).getMaxx();
      }
    }
    LP_ts[1][j].addFirstPiece(new Piece(Track(), Interval(mini, maxi), Cost()));

    for(unsigned int i = 2; i < (n + 1); i++)
    {
      LP_ts[i][j].addFirstPiece(new Piece(Track(), Interval(mini, maxi), Cost()));
      LP_ts[i][j].setUniquePieceCostToInfinity();
    }
    mini = inter.geta();
    maxi = inter.getb();
  }

  ///START STATE CONSTRAINT + add FirstPoint to LP_ts[1]
  ///START STATE CONSTRAINT + add FirstPoint to LP_ts[1]
  std::vector<unsigned int> startState = m_graph.getStartState();
  if(startState.size() != 0)
  {
    for(unsigned int j = 0; j < p; j++)
    {
      if(std::find(startState.begin(), startState.end(), j) == startState.end())
        {LP_ts[1][j].setUniquePieceCostToInfinity();}
      else{LP_ts[1][j].initializeHeadWithFirstPoint(firstData);}
    }
  }
  else
  {
    for(unsigned int j = 0; j < p; j++)
      {LP_ts[1][j].initializeHeadWithFirstPoint(firstData);}
  }
}

//####### gfpop BEGIN #######// //####### gfpop BEGIN #######// //####### gfpop BEGIN #######//
//####### gfpop BEGIN #######// //####### gfpop BEGIN #######// //####### gfpop BEGIN #######//
//####### gfpop BEGIN #######// //####### gfpop BEGIN #######// //####### gfpop BEGIN #######//
//####### gfpop BEGIN #######// //####### gfpop BEGIN #######// //####### gfpop BEGIN #######//

void Omega::gfpop(Data const& data)
{
	Point* myData = data.getVecPt(); // GET the data = vector of Point = myData
  n = data.getn(); // data length
	initialize_LP_ts(myData[0], n); // Initialize LP_ts Piece : size LP_ts (n+1) x p + add first data point

	for(unsigned int t = 1; t < n; t++) // loop for all data point
	{
	  LP_edges_operators(t); // fill_LP_edges. t = newLabel to consider
	  LP_edges_addPointAndPenalty(myData[t]); // Add new data point and penalty
    LP_t_new_multipleMinimization(t); // multiple_minimization
	}
	backtracking();
}

//####### gfpop END #######// //####### gfpop END #######// //####### gfpop END #######//
//####### gfpop END #######// //####### gfpop END #######// //####### gfpop END #######//
//####### gfpop END #######// //####### gfpop END #######// //####### gfpop END #######//
//####### gfpop END #######// //####### gfpop END #######// //####### gfpop END #######//

//####### gfpopTestMode BEGIN #######// //####### gfpopTestMode BEGIN #######// //####### gfpopTestMode BEGIN #######//
//####### gfpopTestMode BEGIN #######// //####### gfpopTestMode BEGIN #######// //####### gfpopTestMode BEGIN #######//

void Omega::gfpopTestMode(Data const& data)
{

  Point* myData = data.getVecPt(); // GET the data = vector of Point = myData
  n = data.getn(); // data length
  initialize_LP_ts(myData[0], n); // Initialize LP_ts Piece : size LP_ts (n+1) x p + add first data point

  for(unsigned int i = 0; i < p; i++)
  {
    //std::cout << "position "<< 1 << " ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"<< std::endl;
    //std::cout << "state "<< i << std::endl;
    LP_ts[1][i].show();
    //std::cout << " TEST ";
    LP_ts[1][i].test();
  }

  for(unsigned int t = 1; t < n; t++) // loop for all data point
  {
    // /*
    //std::cout << std::endl;
    //std::cout << t << "  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
    // */
    LP_edges_operators(t); // fill_LP_edges. t = newLabel to consider
    LP_edges_addPointAndPenalty(myData[t]); // Add new data point and penalty

    // /*
    //std::cout << "  LP_edgesLP_edgesLP_edgesLP_edgesLP_edgesLP_edges "<< t<< std::endl;
    for(unsigned int i = 0; i < q; i++) /// loop for all q edges
    {
      //std::cout << i << "  type " << m_graph.getEdge(i).getConstraint() << "  states " << m_graph.getEdge(i).getState1() << " and " << m_graph.getEdge(i).getState2() << std::endl;
      LP_edges[i].show();
      //std::cout << " TEST ";
      LP_edges[i].test();
    }
    //std::cout << "----------------------------------------------------------------------------------------------------------------------------------"<< std::endl;
    // */
    LP_t_new_multipleMinimization(t); // multiple_minimization

    // /*
    for(unsigned int i = 0; i < p; i++)
    {
      //std::cout << "position "<< t + 1 << " ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"<< std::endl;
      //std::cout << "state "<< i << std::endl;
      LP_ts[t+1][i].show();
      //std::cout << t << " state "<< i << " ";
      LP_ts[t+1][i].test();
    }
    Rcpp::checkUserInterrupt();
    // */

  }
  backtracking();
}

//####### gfpopTestMode END #######// //####### gfpopTestMode END #######// //####### gfpopTestMode END #######//
//####### gfpopTestMode END #######// //####### gfpopTestMode END #######// //####### gfpopTestMode END #######//


//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///
//##### LP_edges_operators #####//////##### LP_edges_operators #####//////##### LP_edges_operators #####///

void Omega::LP_edges_operators(unsigned int t)
{
  for(unsigned int i = 0 ; i < q ; i++) /// loop for all q edges
  {
    // COMMENT: i-th edge = m_graph.getEdge(i)
    // COMMENT: starting state = m_graph.getEdge(i).getState1()
    // COMMENT: t is the label to associate to the constraint
    LP_edges[i].LP_edges_constraint(LP_ts[t][m_graph.getEdge(i).getState1()], m_graph.getEdge(i), t);
  }
}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void Omega::LP_edges_addPointAndPenalty(Point const& pt)
{
  for(unsigned char i = 0; i < q; i++) /// loop for all q edges
  {
    // COMMENT: LP_edges[i] = i-th edge = m_graph.getEdge(i) BECAUSE we need K, a and penalty
    LP_edges[i].LP_edges_addPointAndPenalty(m_graph.getEdge(i), pt);
  }
}

//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///
//##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####//////##### multipleMinimization_LP_edges #####///

void Omega::LP_t_new_multipleMinimization(unsigned int t)
{
  // COMMENT: m_graph was rearranged with increasing integer state2 AND increasing beta penalty
  // COMMENT: LP_ts[t + 1][j] initialized in initialize_LP_ts by addFirstPiece(new Piece(Track(), Interval(mini, maxi), +INFINITY))
  unsigned int k = 0;
  for(unsigned int j = 0 ; j < p; j++)
  {
    while((k < q) && (m_graph.getEdge(k).getState2() == j))
    {
      LP_ts[t + 1][j].LP_ts_Minimization(LP_edges[k]);
      k = k + 1;
      /*
      std::cout << " SHOW ";
      LP_ts[t + 1][j].show();
      std::cout << " TEST ";
      LP_ts[t + 1][j].test();
       */
    }
  }
}


//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///
//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///

void Omega::backtracking()
{
  std::vector< int > changepoints1;
  std::vector< double > parameters1;
  std::vector< int > states1;
  std::vector< int > forced1;

  Interval constrainedInterval; // Interval to fit the constraints

  double* malsp = new double[5];
  double* malsp_temp = new double[5];
  //Interval* nodeConstr = m_graph.nodeConstraints();

  LP_ts[n][0].get_min_argmin_label_state_position_ListPiece(malsp);


  ///////////////////
  /// FINAL STATE ///
  ///////////////////
  unsigned int CurrentState = 0; // Current state
  unsigned int CurrentChgpt = n; // data(1)....data(n). Last data index in each segment
  double CurrentGlobalCost;
  std::vector<unsigned int> endState = m_graph.getEndState();

  // IF no endState, all the states are endstates => we select the best one
  if(endState.size() == 0)
  {
    for (unsigned int j = 1 ; j < p ; j++) // for all p states
    {
      LP_ts[n][j].get_min_argmin_label_state_position_ListPiece(malsp_temp);
      if(malsp_temp[0] < malsp[0]){CurrentState = j; malsp[0] = malsp_temp[0];}
    }
  }
  else ///=> multiple output //to be done
  {
    malsp[0] = INFINITY;
    for (unsigned int j = 0 ; j < endState.size() ; j++) // for all endState available
    {
      LP_ts[n][endState[j]].get_min_argmin_label_state_position_ListPiece(malsp_temp);
      if(malsp_temp[0] < malsp[0]){CurrentState = endState[j]; malsp[0] = malsp_temp[0];}
    }
  }

  ///// with the best state
  LP_ts[n][CurrentState].get_min_argmin_label_state_position_ListPiece(malsp);
  CurrentGlobalCost = malsp[0];
  parameters1.push_back(malsp[1]); // = argmin
  changepoints1.push_back(CurrentChgpt); // = n
  states1.push_back(CurrentState); // = the best state

  /// BACKTRACK
  ///////////////////////////////
  /// previous to FINAL STATE ///
  ///////////////////////////////

  bool out;
  bool boolForced;
  double decay = 0;
  double correction = 1;

  while(malsp[2] > 0) ///while Label > 0
  {
    out = false;
    boolForced = false;
    decay = m_graph.recursiveState(CurrentState);
    if(decay != 1){correction = std::pow(1.0*decay, 1.0*(parameters1.back() - malsp[2] + 1));}else{correction = 1;}

    constrainedInterval = m_graph.buildInterval(malsp[1]*correction, malsp[3], CurrentState, out); ///update out

    CurrentState = malsp[3];
    CurrentChgpt = malsp[2];

    //TO UPDATE: malsp[4] = position
    LP_ts[(int) malsp[2]][(int) malsp[3]].get_min_argmin_label_state_position_onePiece(malsp, (int) malsp[4], constrainedInterval, out, boolForced); ///update boolForced

    //update CurrentGlobalCost and boolForced if argmin on a bound
    CurrentGlobalCost = CurrentGlobalCost - m_graph.findBeta(malsp[3], CurrentState);
    //if(malsp[1] == nodeConstr[CurrentState].geta() || malsp[1] == nodeConstr[CurrentState].getb()){boolForced = true;}

    parameters1.push_back(malsp[1]);
    changepoints1.push_back(CurrentChgpt);
    states1.push_back(CurrentState);
    forced1.push_back(boolForced);
  }

  globalCost.push_back(CurrentGlobalCost);

  ////
  parameters.push_back(parameters1);
  changepoints.push_back(changepoints1);
  states.push_back(states1);
  forced.push_back(forced1);

  delete[] malsp;
  delete[] malsp_temp;
  //delete(nodeConstr);
}


///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###
///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###

void Omega::show()
{
  for(unsigned int i = 0; i < p; i++) /// loop for all q edges
  {
    //std::cout << "  type " << m_graph.getEdge(i).getConstraint() << std::endl;
    //std::cout << "  states " << m_graph.getEdge(i).getState1() << " and " << m_graph.getEdge(i).getState2() << std::endl;
    LP_edges[i].show();
  }

}

