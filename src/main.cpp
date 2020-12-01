#include <iostream>
#include <string>
#include<math.h>

#include<Rcpp.h>

#include"Omega.h"
#include"Cost.h"
#include"ExternFunctions.h"

#include"Data.h"
#include"Graph.h"
#include"Edge.h"

using namespace Rcpp;

// [[Rcpp::export]]
List gfpopTransfer(NumericVector vectData, DataFrame mygraph, std::string type, NumericVector vectWeight, bool testMode)
{
  ///////////////////////////////////////////
  /////////// DATA TRANSFORMATION ///////////
  ///////////////////////////////////////////
  double epsilon = std::pow(10.0,-12.0);

  if(type == "variance")
  {
    double mean = 0;
    for(int i = 0; i < vectData.size(); i++){mean = mean + vectData[i];}
    mean = mean/vectData.size();
    for(int i = 0; i < vectData.size(); i++){vectData[i] = vectData[i] - mean; if(vectData[i] == 0){vectData[i] = epsilon;}}
  }

  if(type == "poisson")
  {
    for(int i = 0; i < vectData.size(); i++){if(vectData[i] < 0){throw std::range_error("There are some negative data");}}
    //for(int i = 0; i < vectData.size(); i++){if(vectData[i]  > floor(vectData[i])){throw std::range_error("There are some non-integer data");}}
  }

  if(type == "exp")
  {
    for(int i = 0; i < vectData.size(); i++){if(vectData[i] <= 0){throw std::range_error("Data has to be all positive");}}
  }

  if(type == "negbin")
  {
    unsigned int windowSize = 100;
    unsigned int k = vectData.size() / windowSize;
    double mean = 0;
    double variance = 0;
    double disp = 0;

    for(unsigned int j = 0; j < k; j++)
    {
      mean = 0;
      variance = 0;
      for(unsigned int i = j * windowSize; i < (j + 1)*windowSize; i++){mean = mean + vectData[i];}
      mean = mean/windowSize;
      for(unsigned int i =  j * windowSize; i < (j + 1)*windowSize; i++){variance = variance + (vectData[i] - mean) * (vectData[i] - mean);}
      variance = variance/(windowSize - 1);
      disp = disp  + (mean * mean / (variance - mean));
    }
    disp = disp/k;
    for(int i = 0; i < vectData.size(); i++){vectData[i] = vectData[i]/disp; if(vectData[i] == 0){vectData[i] = epsilon/(1- epsilon);}}
  }

  // BEGIN TRANSFERT into C++ objects  // BEGIN TRANSFERT into C++ objects  // BEGIN TRANSFERT into C++ objects
  // BEGIN TRANSFERT into C++ objects  // BEGIN TRANSFERT into C++ objects  // BEGIN TRANSFERT into C++ objects
  // DATA AND GRAPH

  /////////////////////////////////
  /////////// DATA COPY ///////////
  /////////////////////////////////
  Data data = Data();
  data.copy(vectData, vectWeight, vectData.length(), vectWeight.length());

  //////////////////////////////////
  /////////// GRAPH COPY ///////////
  //////////////////////////////////

  Graph graph = Graph();
  Edge newedge;

  ///9 variables in mygraph
  Rcpp::IntegerVector state1 = mygraph["state1"];
  Rcpp::IntegerVector state2 = mygraph["state2"];
  Rcpp::CharacterVector typeEdge = mygraph["type"];
  Rcpp::NumericVector parameter = mygraph["parameter"];
  Rcpp::NumericVector penalty = mygraph["penalty"];
  Rcpp::NumericVector KK = mygraph["K"];
  Rcpp::NumericVector aa = mygraph["a"];
  Rcpp::NumericVector minn = mygraph["min"];
  Rcpp::NumericVector maxx = mygraph["max"];

  for(int i = 0 ; i < mygraph.nrow(); i++)
    {graph << Edge(state1[i], state2[i], typeEdge[i], fabs(parameter[i]), penalty[i], fabs(KK[i]), fabs(aa[i]), minn[i], maxx[i]);}

  if(testMode == TRUE){graph.show();} ///TESTMODE

  // END TRANSFERT into C++ objects  // END TRANSFERT into C++ objects  // END TRANSFERT into C++ objects
  // END TRANSFERT into C++ objects  // END TRANSFERT into C++ objects  // END TRANSFERT into C++ objects

  /////////////////////////////////////////////
  /////////// COST FUNCTION LOADING ///////////
  /////////////////////////////////////////////

  cost_coeff = coeff_factory(type);
  cost_eval = eval_factory(type);

  cost_min = min_factory(type);
  cost_minInterval = minInterval_factory(type);
  cost_argmin = argmin_factory(type);
  cost_argminInterval = argminInterval_factory(type);
  cost_argminBacktrack = argminBacktrack_factory(type);

  cost_shift = shift_factory(type);
  cost_interShift = interShift_factory(type);
  cost_expDecay = expDecay_factory(type);
  cost_interExpDecay = interExpDecay_factory(type);

  cost_intervalInterRoots = intervalInterRoots_factory(type);
  cost_age = age_factory(type);
  cost_interval = interval_factory(type);

  /////////////////////////////
  /////////// OMEGA ///////////
  /////////////////////////////

  Omega omega(graph);
  if(testMode == FALSE){omega.gfpop(data);}else{omega.gfpopTestMode(data);}

  /////////////////////////////
  /////////// RETURN //////////
  /////////////////////////////

  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["states"] = omega.GetStates(),
    _["forced"] = omega.GetForced(),
    _["param"] = omega.GetParameters(),
    _["cost"] = omega.GetGlobalCost()
);

  return res;
}
