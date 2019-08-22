//#include <iostream>
#include <string>
#include<math.h>

#include<Rcpp.h>

#include"Omega.h"

#include"Data.h"
#include"Graph.h"
#include"Edge.h"
#include"Bound.h"
#include"Robust.h"

using namespace Rcpp;

// [[Rcpp::export]]
List gfpopTransfer(NumericVector vectData, DataFrame mygraph, std::string type, NumericVector vectWeight)
{
  // BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT
  // BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT// BEGIN TRANSFERT
  // BEGIN TRANSFERT INTO C++ OBJETS

  ///////////
  /////////// DATA LOADING
  ///////////
  Data data = Data(); // in any case, add a file name. by default = "nofile"
  int n = vectData.length();
  int nw = vectWeight.length();
  data.copy(vectData, vectWeight, n, nw);
  //data.show(false);

  ///////////
  /////////// GRAPH
  ///////////

  Graph graph = Graph();
  Edge newedge;

  Rcpp::IntegerVector state1 = mygraph["state1"];
  Rcpp::IntegerVector state2 = mygraph["state2"];
  Rcpp::CharacterVector typeEdge = mygraph["type"];
  Rcpp::NumericVector penalty = mygraph["penalty"];
  Rcpp::NumericVector parameter = mygraph["parameter"];

  for (int i = 0 ; i < mygraph.nrow(); i++)
    {graph << Edge(penalty[i], state1[i], state2[i], typeEdge[i], parameter[i]);}
  ///Include start and end states if specified

  std::string graphType = graph.getType();
  //Rcout << "graphType : " << graphType << std::endl;
  graph.show();
  if(graph.AreVerticesCompatible() == false){Rcout << "The vertices must be labeled by integers from 0 to S (an integer)" << std::endl; return(0);}


  ///////////TO DELETE
  Bound bound = Bound(0, 0, false);
  Robust robust = Robust(1000,1000);

  // END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT
  // END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT// END TRANSFERT

  ///////////
  /////////// OMEGA
  ///////////

  Omega omega(graph, bound, robust);

  /// FPOP 1D GRAPH ALGORITHM. Different call with respect to the complexity of the graph.
  if(graphType == "std"){omega.fpop1d_graph_std(data);}
  if(graphType == "isotonic"){omega.fpop1d_graph_isotonic(data);}
  if(graphType == "pava"){omega.pava(data);}
  if(graphType == "complex"){omega.fpop1d_graph_complex(data);}

  ///////////
  /////////// RETURN
  ///////////

  List res = List::create(
    _["changepoints"] = omega.GetChangepoints(),
    _["states"] = omega.GetStates(),
    _["forced"] = omega.GetForced(),
    _["param"] = omega.GetParameters(),
    _["cost"] = omega.GetGlobalCost()
  );

  return res;
}


