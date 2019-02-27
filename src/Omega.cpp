#include "Omega.h"

#include "CostGauss.h"

#include "Piece.h"

#include<iostream>
#include <iomanip> ///Set Precision scientific

#include<vector>
#include<typeinfo>
#include <stdlib.h>

#include <fstream> ///write in a file
#include<sstream> ///for conversion int to string : function ostringstream

#include <algorithm>

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

Omega::Omega(Graph graph, Bound bound, Robust robust) : m_graph(graph), m_bound(bound), m_robust(robust)
{
	p = graph.nb_states();
	q = graph.nb_edges();
  double mini = m_bound.getm();
  double maxi = m_bound.getM();

	Piece **initRow = new Piece*[p];
  for(unsigned char i = 0 ; i < p ; i++){initRow[i] = new Piece(Track(), Interval(mini, maxi), CostGauss());}
	Q_ts.push_back(initRow); ///push into Q_ts as the first element (first element of the vector Q_ts)

	Q_edges = new Piece*[q];
	for(unsigned char i = 0 ; i < q ; i++){Q_edges[i] = new Piece(Track(), Interval(mini, maxi), CostGauss());} ///because we use delete Q_edges in loops
	Q_s_temp = new Piece*[p];
  for(unsigned char i = 0 ; i < p ; i++){Q_s_temp[i] = new Piece(Track(), Interval(mini, maxi), CostGauss());} ///because we use delete Q_s_temp in loops
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Omega::~Omega()
{
	for(int i = 0; i < Q_ts.size(); i++){for(unsigned char j = 0; j < p; j++){delete(Q_ts[i][j]);}delete [] Q_ts[i]; Q_ts[i] = NULL;}

  for(unsigned char i = 0; i < q; i++){if(Q_edges[i] != NULL){delete(Q_edges[i]);}}
  delete [] Q_edges;
  Q_edges = NULL;

  for(unsigned char i = 0; i < p; i++){delete(Q_s_temp[i]);}
  delete [] Q_s_temp;
  Q_s_temp = NULL;
}

//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//

std::vector< int > Omega::GetChangepoints() const{return(changepoints);}
std::vector< double > Omega::GetMeans() const{return(means);}
std::vector< int > Omega::GetStates() const{return(states);}
std::vector< int > Omega::GetForced() const{return(forced);}

int Omega::GetN() const{return(Q_ts.size());}
double Omega::GetGlobalCost() const{return(globalCost);}

void Omega::test(Data const& data)
{
  Q_edges = new Piece*[data.getn()];
  for(unsigned char i = 0 ; i < p ; i++){Q_edges[i] = new Piece(Track(), Interval());}
}


//####### fpop1d_graph #######// //####### fpop1d_graph #######// //####### fpop1d_graph #######//
//####### fpop1d_graph #######// //####### fpop1d_graph #######// //####### fpop1d_graph #######//
//####### fpop1d_graph #######// //####### fpop1d_graph #######// //####### fpop1d_graph #######//
//####### fpop1d_graph #######// //####### fpop1d_graph #######// //####### fpop1d_graph #######//

//####### complex #######////####### complex #######////####### complex #######//
//####### complex #######////####### complex #######////####### complex #######//

void Omega::fpop1d_graph_complex(Data const& data)
{
	Point* myData = data.getVecPt(); ///GET the data/// get the vector of Points = myData

	///
	/// Initial functional cost (with startState constraint)
	///
	///add a new row in Q_ts = copy the previous one
	Piece** Qt_new = new Piece*[p];
	for (unsigned char i = 0; i < p ; i++){Qt_new[i] = Q_ts.back()[i] -> copy();}
	Q_ts.push_back(Qt_new);

  addPointQ_t(myData[0]);
  std::vector<int> startState = m_graph.getStartState();
  if(startState.size() != 0){for(int i = 0; i < p; i++){if(std::find(startState.begin(), startState.end(), i) == startState.end()){Q_ts.back()[i] -> addConstant(INFINITY);}}}

  for(unsigned int t = 1; t < data.getn(); t++) /// loop for all edges
  {
    fillQ_edges(t); ///fillQ_edges. t = newLabel to consider
    copyQt(); ///add a new row in Q_ts = copy the previous one
    multiple_minimization(); ///multiple_minimization
    addPointQ_t(myData[t]); ///Add new data point
  }
	backtracking();
}


//####### isotonic #######////####### isotonic #######////####### isotonic #######//
//####### isotonic #######////####### isotonic #######////####### isotonic #######//

void Omega::fpop1d_graph_isotonic(Data const& data)
{
  Point* myData = data.getVecPt(); ///GET the data///get the vector of Points = myData

  ///
  /// new Pointers
  ///
  Piece* Q_up = new Piece(Q_edges[0] -> getTrack(), Q_edges[0] -> getInterval(), CostGauss());
  Piece* Q_t0 = new Piece(Q_edges[0] -> getTrack(), Q_edges[0] -> getInterval(), CostGauss());
  std::vector<Piece*> Q_t;
  Q_t.push_back(Q_t0);

  ///
  /// BUILD the array currentMin
  ///
  double* currentMin = new double[data.getn()];
  double current = myData[data.getn() - 1].y; ///= last data
  if(current > m_bound.getM()){current = m_bound.getM();}
  currentMin[data.getn() - 1] = current;

  for(int t = data.getn() - 2 ; t > - 1; t--)
  {
    if(myData[t].y < current){current = myData[t].y;}
    if(current > m_bound.getM()){current = m_bound.getM();}
    currentMin[t] = current;
  }
  double newLeftBound = currentMin[0]; ///the minimal value in myData -> the value for the new left bound at each iteration


  ///
  /// ISOTONIC parameters
  ///
  double betaUp = m_graph.getEdge(1).getBeta();
  double parameterUp = m_graph.getEdge(1).getParameter();

  ///
  /// add point in recursion
  ///
  for(unsigned int t = 0 ; t < data.getn(); t++)
  {
		newLeftBound = Q_t.back() -> newBound(currentMin[t]); /// left bound update = update partialArgmin (newLeftBound <= currentMin[t])
    Q_t.push_back(Q_t.back() -> copyIsotonic(newLeftBound)); ///new element in Q_t => copy Q_t with the updated left bound

    delete(Q_up); /// reset the vector of Piece Q_edges[0]
    Q_up = NULL;
    Q_up = Q_t.back() -> operator_down(t, -1); /// parentStateLabel = -1 (a unique state => "no state to consider")
    Q_up -> addConstant(betaUp); ///add beta to Q_up
    if(parameterUp > 0){Q_up = Q_up -> shift_right(parameterUp, m_bound.getM());} ///if parameter > 0 => shift right

    Q_t.back() = Q_t.back() -> min_function(Q_up, m_bound.getM()); /// minimum opertor : Q_up VS Q_t
    Q_t.back() -> addPoint(myData[t], m_robust);

    /*std::ostringstream ttemp;  //temp as in temporary
    ttemp << t;
    std::string fileQup = "Rtxt/Qup_" + ttemp.str();
    std::string fileQt = "Rtxt/Qt_" + ttemp.str();

    fileQup = fileQup + ".txt";
    fileQt = fileQt + ".txt";

    std::ofstream fichier1(fileQup.c_str(),std::ios::out | std::ios::trunc);
    fichier1.precision(9);
    Piece* tmp1 = Q_up -> copy();
    fichier1 >> tmp1;
    fichier1.close();
    std::ofstream fichier2(fileQt.c_str(),std::ios::out | std::ios::trunc);
    fichier2.precision(9);
    Piece* tmp2 = Q_t.back() -> copy();
    fichier2 >> tmp2;
    fichier2.close();*/

  }
  backtrackingIsotonic(Q_t);

  for(unsigned int i = 0 ; i < Q_t.size(); i++){delete(Q_t[i]);}
  delete [] currentMin;
  delete(Q_up);
}



//####### std #######////####### std #######////####### std #######//
//####### std #######////####### std #######////####### std #######//

void Omega::fpop1d_graph_std(Data const& data)
{
  double beta = m_graph.getEdge(1).getBeta();

  Point* myData = data.getVecPt(); ///GET the data///get the vector of Points = myData

  int lastlabel;
  double argmini;
  std::vector< int > temp_changepoints; ///vector of changepoints build by fpop
  std::vector< double > temp_means; ///vector of means build by fpop

  ///2 cases.
  ///bound constrainted = true
  if(m_bound.getIsConstrained() == true)
  {
    for(unsigned int t = 0 ; t < data.getn(); t++) /// loop for all edges
    {
      Q_s_temp[0] -> addPoint(myData[t], m_robust);
      delete(Q_edges[0]);
      Q_edges[0] = NULL;
      Q_edges[0] = Q_s_temp[0] -> operator_stdConstr_min_argmin(t, lastlabel, argmini, m_bound); /// = one Piece constant at global min.

      temp_changepoints.push_back(lastlabel);
      temp_means.push_back(argmini);

      Q_edges[0] ->   getRefCost() += beta; /// add beta (the penalty)
      Q_s_temp[0] = Q_s_temp[0] -> min_function(Q_edges[0], m_bound.getM());
    }
  }

  ///bound constrainted = false
  if(m_bound.getIsConstrained() == false)
  {
    for(unsigned int t = 0 ; t < data.getn(); t++) /// loop for all edges
    {
      Q_s_temp[0] -> addPoint(myData[t], m_robust);
      delete(Q_edges[0]);
      Q_edges[0] = NULL;
      Q_edges[0] = Q_s_temp[0] -> operator_std_min_argmin(t, lastlabel, argmini, m_bound); /// = one Piece constant at global min.

      temp_changepoints.push_back(lastlabel);
      temp_means.push_back(argmini);

      Q_edges[0] -> getRefCost() += beta; /// add beta (the penalty)
      Q_s_temp[0] = Q_s_temp[0] -> min_function(Q_edges[0], m_bound.getM());
    }
  }

  /// BACKTRACK
  int length = temp_changepoints.size() - 1;
  int position = temp_changepoints[length];

  ///Find globalCost
  std::vector<double> response = Q_s_temp[0] -> get_min_argmin_label_state_position_final();

  ///Fill Omega VARIABLES
  globalCost = response[0];
  changepoints.push_back(data.getn());

  changepoints.push_back(position + 1);
  means.push_back(temp_means[length]);

  while(position > 0)
  {
    means.push_back(temp_means[position]);
    position = temp_changepoints[position];
    if(position != 0){changepoints.push_back(position + 1);}
  }

}




//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///
//##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####/// ///##### SUBFUNCTIONS #####///

//##### copyQt #####//////##### copyQt #####//////##### copyQt #####///
//##### copyQt #####//////##### copyQt #####//////##### copyQt #####///

void Omega::copyQt()
{
	//Piece** Qt = Q_ts.back();
	Piece** Qt_new = new Piece*[p];
	//for (unsigned char i = 0; i < p ; i++){Qt_new[i] = Qt[i] -> copy();}
	Q_ts.push_back(Qt_new);
}


//##### fillQ_edges #####//////##### fillQ_edges #####//////##### fillQ_edges #####///
//##### fillQ_edges #####//////##### fillQ_edges #####//////##### fillQ_edges #####///

void Omega::fillQ_edges(int newLabel)
{
	int s1;
	for (unsigned int i = 0 ; i < q ; i++) /// loop for all edges
	{
    delete(Q_edges[i]);///DELETE Q_edges[i]
		Edge edge = m_graph.getEdge(i);
    s1 = edge.getState1();   /// starting state
    Q_edges[i] = Q_ts.back()[s1] -> edge_constraint(edge, newLabel, m_bound);
	}
}

//##### multiple_minimization #####//////##### multiple_minimization #####//////##### multiple_minimization #####///
//##### multiple_minimization #####//////##### multiple_minimization #####//////##### multiple_minimization #####///

void Omega::multiple_minimization()
{
  int j = 0;
  ///Q_s_temp vs Q_ts minimization
  for (unsigned int i = 0 ; i < p; i++)
  {
    //delete(Q_ts.back()[i]);
    Q_ts.back()[i] = Q_edges[j] -> copy();

    while((j + 1 < q) && (m_graph.getEdge(j + 1).getState2() == i))
    {
      //Piece* Q_edgesCopy = Q_edges[j + 1] -> copy();
      Q_ts.back()[i] = Q_ts.back()[i] -> min_function(Q_edges[j + 1], m_bound.getM()); ///
      //delete(Q_edgesCopy);
      j = j + 1;
    }
    j = j + 1;
  }
}

//##### addPointQ_t #####//////##### addPointQ_t #####//////##### addPointQ_t #####///
//##### addPointQ_t #####//////##### addPointQ_t #####//////##### addPointQ_t #####///

void Omega::addPointQ_t(Point pt)
{
	for (unsigned char i = 0; i < p; i++){Q_ts.back()[i] -> addPoint(pt, m_robust);}
}


//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///
//##### backtracking #####//////##### backtracking #####//////##### backtracking #####///

void Omega::backtracking()
{
  Interval constrainedInterval; ///Interval to fit the constraints
  ///
  /// malsp = Min_Argmin_Label_State_Position_Final
  ///
  std::vector<double> malsp = Q_ts.back()[0] -> get_min_argmin_label_state_position_final();
  std::vector<double> malsp_temp;

  ///
  ///FINAL STATE
  ///
  int CurrentState = 0; ///Current state
  int CurrentChgpt = Q_ts.size() - 1; /// data(1)....data(n). Last data index in each segment
  std::vector<int> endState = m_graph.getEndState();

  if(endState.size() == 0)
  {
    for (int j = 1 ; j < p ; j++) ///for all states
    {
      malsp_temp = Q_ts.back()[j] -> get_min_argmin_label_state_position_final();
      if(malsp_temp[0] < malsp[0]){CurrentState = j; malsp[0] = malsp_temp[0];}
    }
  }
  else
  {
    for (int j = 0 ; j < endState.size() ; j++) ///for all states
    {
      malsp_temp = Q_ts.back()[endState[j]] -> get_min_argmin_label_state_position_final();
      if(malsp_temp[0] < malsp[0]){CurrentState = endState[j]; malsp[0] = malsp_temp[0];}
    }
  }

  malsp = Q_ts.back()[CurrentState] -> get_min_argmin_label_state_position_final();
  globalCost = malsp[0];

  //std::cout << std::endl;
  //std::cout << "BACKTRACKING" << std::endl;
  //std::cout<< "CurrentChgpt : "<< CurrentChgpt << " -- " << "CurrentState : "<< CurrentState << std::endl;
  //Q_ts.back()[CurrentState] -> show(); ///Q_ts LABEL,STATE
  ///
  ///SHOW AND SAVE DATA
  ///
  //std::cout << "Analysis" << std::endl;
  //std::cout<< "     lab : "<< malsp[2] << std::endl;
  //std::cout<< "     state : "<< malsp[3] << std::endl;
  //std::cout<< "     position : " << malsp[4] << std::endl;
  //std::cout<< "     min : " << malsp[0] << std::endl;
  //std::cout<< "     arg : "<< malsp[1] << std::endl;

  means.push_back(malsp[1]);
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
    if(decay != 1){correction = std::pow(decay, means.back() - malsp[2] + 1);}

    constrainedInterval = m_graph.buildInterval(malsp[1]*correction, malsp[3], CurrentState, out); ///update out
    CurrentState = malsp[3];
    CurrentChgpt = malsp[2];

    malsp = Q_ts[malsp[2]][(int) malsp[3]] -> get_min_argmin_label_state_position((int) malsp[4], constrainedInterval, out, boolForced, m_bound.getIsConstrained()); ///update boolForced

    if(malsp[1] > m_bound.getM()){malsp[1] = m_bound.getM(); boolForced = true;}
    if(malsp[1] < m_bound.getm()){malsp[1] = m_bound.getm(); boolForced = true;}

    //std::cout << std::endl;
    //std::cout<< "CurrentChgpt : " << CurrentChgpt << " -- " << "CurrentState : " << CurrentState << std::endl;
    //Q_ts[CurrentChgpt][CurrentState] -> show();
    //constrainedInterval.show();
    ///
    ///SHOW AND SAVE DATA
    ///
    //std::cout << "Analysis" << std::endl;
    //std::cout<< "     lab : " << malsp[2] << std::endl;
    //std::cout<< "     state : " << malsp[3] << std::endl;
    //std::cout<< "     position : " << malsp[4] << std::endl;
    //std::cout<< "     min : " << malsp[0] << std::endl;
    //std::cout<< "     arg : " << malsp[1] << std::endl;

    means.push_back(malsp[1]);
    changepoints.push_back(CurrentChgpt);
    states.push_back(CurrentState);
    forced.push_back(boolForced);
  }
  //std::cout<<"DONE"<<std::endl;
}




//##### backtrackingIsotonic #####//////##### backtrackingIsotonic #####//////##### backtrackingIsotonic #####///
//##### backtrackingIsotonic #####//////##### backtrackingIsotonic #####//////##### backtrackingIsotonic #####///

void Omega::backtrackingIsotonic(std::vector<Piece*> const& Q_t)
{
  ///
  /// malsp = Min_Argmin_Label_State_Position_Final
  ///

  std::vector<double> malsp = Q_t.back() -> get_min_argmin_label_state_position_final();
  globalCost = malsp[0];
  //std::cout<< "     min : " << malsp[0] << std::endl;

  int CurrentChgpt = Q_t.size() - 1; /// data(1)....data(n). Last data index in each segment

  //std::cout << std::endl;
  //std::cout << "BACKTRACKING" << std::endl;
  //std::cout<< "CurrentChgpt : "<< CurrentChgpt << std::endl;
  //Q_t.back() -> show(); ///Q_ts LABEL,STATE
  ///
  ///SHOW AND SAVE DATA
  ///
  //std::cout << "Analysis" << std::endl;
  //std::cout<< "     lab : "<< malsp[2] << std::endl;
  //std::cout<< "     min : " << malsp[0] << std::endl;
  //std::cout<< "     arg : "<< malsp[1] << std::endl;

  means.push_back(malsp[1]);
  changepoints.push_back(CurrentChgpt);
  states.push_back(0); ///the only state is vertex state 0


  /// BACKTRACK
  ///
  ///BEFORE FINAL STATE
  ///
  bool boolForced = false;

  while( malsp[2] > 0)  ///while Label > 0
  {
    ///
    ///BACKTRACK
    ///
    boolForced = false;
    CurrentChgpt = malsp[2];
    malsp = Q_t[malsp[2]] -> get_min_argmin_label(malsp[1] - m_graph.getEdge(1).getParameter(), boolForced, m_bound.getIsConstrained());

    //malsp = Q_t[malsp[2]] -> get_min_argmin_label(INFINITY, boolForced, m_bound.getIsConstrained());
    if(malsp[1] > m_bound.getM()){malsp[1] = m_bound.getM(); boolForced = true;}
    if(malsp[1] < m_bound.getm()){malsp[1] = m_bound.getm(); boolForced = true;}

    //std::cout <<  malsp[1] << std::endl;
    //std::cout<< "     min : " << malsp[0] << std::endl;

    //std::cout<< "CurrentChgpt : " << CurrentChgpt << std::endl;
    //Q_t[CurrentChgpt] -> show();
    ///
    ///SHOW AND SAVE DATA
    ///
    //std::cout << "Analysis" << std::endl;
    //std::cout<< "     lab : " << malsp[2] << std::endl;
    //std::cout<< "     min : " << malsp[0] << std::endl;
    //std::cout<< "     arg : " << malsp[1] << std::endl;

    means.push_back(malsp[1]);
    changepoints.push_back(CurrentChgpt);
    states.push_back(0);
    forced.push_back(boolForced);
  }

  //std::cout<<"DONE"<<std::endl;
  //for(unsigned int i = 0 ; i < Q_t.size(); i++){delete(Q_t[i]);}

}


//##### saveInFiles #####/// ///##### saveInFiles #####/// ///##### saveInFiles #####/// ///##### saveInFiles #####///
//##### saveInFiles #####/// ///##### saveInFiles #####/// ///##### saveInFiles #####/// ///##### saveInFiles #####///


//##### save_fillQ_edges #####//////##### save_fillQ_edges #####//////##### save_fillQ_edges #####//////##### save_fillQ_edges #####///
//##### save_fillQ_edges #####//////##### save_fillQ_edges #####//////##### save_fillQ_edges #####//////##### save_fillQ_edges #####///

void Omega::save_Q_ts_Q_edges(int t) const
{
  std::ostringstream ttemp;  //temp as in temporary
	ttemp << t;
  std::string fileQ_ts = "Rtxt/Q_ts_" + ttemp.str();
	std::string fileQ_edges = "Rtxt/Q_edges_" + ttemp.str();

	fileQ_ts = fileQ_ts + "_";
  fileQ_edges = fileQ_edges + "_";

  int s1 =  0;

	for (int i = 0; i < q ; i++)
  {
    ///Q_ts

    Edge edge = m_graph.getEdge(i);
    s1 = edge.getState1();   /// starting state

    std::string filenameQ_ts;
    std::ostringstream temp1;  //temp as in temporary
    temp1 << i ;
    filenameQ_ts = fileQ_ts + temp1.str();
    filenameQ_ts += ".txt";
    //std::cout << filenameQ_ts << std::endl;

    std::ofstream fichier1(filenameQ_ts.c_str(),std::ios::out | std::ios::trunc);
    fichier1.precision(9);
    Piece* tmp1 = Q_ts.back()[s1] -> copy();
    fichier1 >> tmp1;
    fichier1.close();

    ///Qedges

    std::string filenameQ_edges;
    std::ostringstream temp2;  //temp as in temporary
    temp2 << i;
    filenameQ_edges = fileQ_edges + temp2.str();
    filenameQ_edges += ".txt";
    //std::cout<< filenameQ_edges <<std::endl;

    std::ofstream fichier2(filenameQ_edges.c_str(),std::ios::out | std::ios::trunc);
    fichier2.precision(9);
    Piece* tmp2 = Q_edges[i] -> copy();
    fichier2 >> tmp2;
    fichier2.close();
  }

}


//##### save_Q_s_temp_Q_ts #####//////##### save_Q_s_temp_Q_ts #####//////##### save_Q_s_temp_Q_ts #####///
//##### save_Q_s_temp_Q_ts #####//////##### save_Q_s_temp_Q_ts #####//////##### save_Q_s_temp_Q_ts #####///


void Omega::save_Q_s_temp_Q_ts(int t) const
{
  std::ostringstream ttemp;  //temp as in temporary
	ttemp << t;
  std::string fileQ_s_temp = "Rtxt/Q_s_temp_" + ttemp.str();
	std::string fileQ_tsNEW = "Rtxt/Q_tsNEW_" + ttemp.str();

	fileQ_s_temp = fileQ_s_temp + "_";
  fileQ_tsNEW = fileQ_tsNEW + "_";

	for (int i = 0; i < p ; i++)
  {
    ///Q_s_temp

    std::string filenameQ_s_temp;
    std::ostringstream temp1;  //temp as in temporary
    temp1 << i ;
    filenameQ_s_temp = fileQ_s_temp + temp1.str();
    filenameQ_s_temp += ".txt";
    //std::cout << filenameQ_s_temp << std::endl;

    std::ofstream fichier1(filenameQ_s_temp.c_str(),std::ios::out | std::ios::trunc);
    fichier1.precision(9);
    Piece* tmp1 = Q_s_temp[i];
    tmp1 -> save(fichier1);
    fichier1.close();

    ///Q_tsNEW

    std::string filenameQ_tsNEW;
    std::ostringstream temp2;  //temp as in temporary
    temp2 << i;
    filenameQ_tsNEW = fileQ_tsNEW + temp2.str();
    filenameQ_tsNEW += ".txt";
    //std::cout<< filenameQ_tsNEW <<std::endl;

    std::ofstream fichier2(filenameQ_tsNEW.c_str(),std::ios::out | std::ios::trunc);

    fichier2.precision(9);
    Piece* tmp2 = Q_ts.back()[i];
    tmp2 -> save(fichier2);

    fichier2.close();
  }
}



///###///###///###///###///###///###///###///###///###///###///###///###///###///###///###

std::ostream &operator<<(std::ostream &s, const Omega &om)
{
  std::vector< int > chpt = om.GetChangepoints();
  std::vector< double > means = om.GetMeans();
  std::vector< int > states = om.GetStates();
  std::vector< int > forced = om.GetForced();
  s << " n : " << om.GetN()-1 << std::endl;
  for (int i = chpt.size()-1; i > -1; i--){s << " ** " << chpt[i];}
  s << std::endl;
  for (int i = means.size()-1; i > -1; i--){s << " ** " << means[i];}
  s << std::endl;
  for (int i = states.size()-1; i > -1; i--){s << " ** " << states[i];}
  s << std::endl;
  for (int i = forced.size()-1; i > -1; i--){s << " ** " << forced[i];}
  s << std::endl;
  s << "globalCost: "<< om.GetGlobalCost() << std::endl;
  return s;
}


