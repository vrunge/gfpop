#include "Piece.h"

#include "termcolor.h"

#include <math.h>
#include <stdlib.h>
#include<iostream>
#include<typeinfo>

#include <fstream> ///write in a file
#include<vector>

Piece::Piece(){m_info = Track(); m_interval = Interval(); m_cost = CostGauss(); nxt = NULL;}
Piece::Piece(Track const& info, Interval const& inter, CostGauss const& cost)
{
  m_info = info;
  m_interval = inter;
  m_cost = cost;
  nxt = NULL;
}

Piece::Piece(const Piece* piece)
{
  m_info = piece -> m_info;
  m_interval = piece -> m_interval;
  m_cost = piece -> m_cost;
  nxt = NULL;
}


//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//

Piece::~Piece()
{
  delete(nxt);
  nxt = NULL;
}


//####### accessors #######////####### accessors #######////####### accessors #######//
//####### accessors #######////####### accessors #######////####### accessors #######//


Track Piece::getTrack()const {return(m_info);}
Interval Piece::getInterval() const {return(m_interval);}
CostGauss Piece::getCost() const {return(m_cost);}

CostGauss& Piece::getRefCost(){return(m_cost);}


//####### global change #######////####### global change #######////####### global change #######//
//####### global change #######////####### global change #######////####### global change #######//


int Piece::length()
{
	Piece* tmp = this;
  int res = 0;
	while(tmp != NULL)
		{
	  res = res + 1;
	  tmp = tmp -> nxt;
		}
	return(res);
}



void Piece::testEmptyInterval()
{
  Piece* tmp = this;
  while(tmp != NULL)
  {
   if(tmp -> m_interval.isEmpty())
     {
        //std::cout << "EMPTY INTERVAL"<< std::endl; exit(0);
     }
    tmp = tmp -> nxt;
  }
}

Piece* Piece::copy()
{
	Piece* tmp = this;
	Piece* firstElement = new Piece(tmp);

	Piece* copyPiece = firstElement;

  tmp = tmp -> nxt;

	while(tmp != NULL)
  {
    copyPiece -> nxt = new Piece(tmp);
    copyPiece = copyPiece -> nxt;
    tmp = tmp -> nxt;
  }

	return(firstElement);
}




Piece* Piece::copy(int& length)
{
	Piece* tmp = this;
	Piece* Qnew = NULL;

	if(tmp != NULL)
		{
      length = length + 1;
			Qnew = new Piece(tmp);
			Qnew -> nxt = tmp -> nxt -> copy(length);
			//delete(tmp);
		}
	return(Qnew);
}




Piece* Piece::copyIsotonic(double newLeftBound)
{
	Piece* tmp = this;
	Piece* Qnew = NULL;

	if(tmp != NULL)
  {
    if(tmp -> m_interval.getb() < newLeftBound)
    {
      Qnew = tmp -> nxt -> copyIsotonic(newLeftBound);
    }
    else
    {
      Qnew = new Piece(tmp);
      Qnew -> nxt = tmp -> nxt -> copyIsotonic(newLeftBound);
      //delete(tmp);
    }
  }
	return(Qnew);
}





Piece* Piece::reverse()
{
  Piece* prev =  NULL;
  Piece* current = this;
  Piece* next = current;

  while(current != NULL)
  {
    current -> m_cost.axisSymmetry();
    current -> m_interval.axisSymmetry();

    next  = current -> nxt;
    current -> nxt = prev;
    prev = current;
    current = next;
  }

  return(prev);
}

Piece* Piece::reverse(int length)
{
  Piece* prev =  NULL;
  Piece* current = this;
  Piece* next = current;

  while(current != NULL)
  {
    current -> m_cost.axisSymmetry();
    current -> m_interval.axisSymmetry();
    current -> m_info.axisSymmetry(length);

    next  = current -> nxt;
    current -> nxt = prev;
    prev = current;
    current = next;
  }

  return(prev);
}


void Piece::opposition()  /// cost <- - cost
{
  Piece* Q_tmp = this;
  while(Q_tmp != NULL)
  {
    Q_tmp -> m_cost.opposition();
    Q_tmp  = Q_tmp -> nxt;
  }
  //delete(Q_tmp);
}




//####### Piece_min_argmin #######////####### Piece_min_argmin #######////####### Piece_min_argmin #######//
//####### Piece_min_argmin #######////####### Piece_min_argmin #######////####### Piece_min_argmin #######//


std::vector<double> Piece::Piece_min_argmin()
{
  std::vector<double> response(2,0); // = {0,0}
  //response[0] = -1;
  //response[1] = -1;

  double argmin = getCost().arg_minimum();

  if(getInterval().isInside(argmin) == true)
    {
      response[1] = argmin;
      response[0] = getCost().minimum();
    }
    else
    {
      if(getInterval().geta() > argmin)
      {
        response[1] = getInterval().geta();
        response[0] =  getCost().point_eval(response[1]);
      }
      else
      {
        response[1] = getInterval().getb();
        response[0] = getCost().point_eval(response[1]);
      }
    }

  return(response);
}



//####### newBound #######// //####### newBound #######// //####### newBound #######//
//####### newBound #######// //####### newBound #######// //####### newBound #######//

double Piece::newBound(double currentDataMin)
{
	Piece* tmp = this;
  double newLeftBound = tmp -> m_interval.geta(); ///newLeftBound = new left bound to create
  double min_temp = tmp -> m_cost.point_eval(newLeftBound);
  std::vector<double> min_argmin;

  while(tmp -> m_interval.getb() < currentDataMin)
  {
    min_argmin = tmp -> Piece_min_argmin(); //min and argmin in Interval
    if(min_argmin[0] < min_temp)
    {
      min_temp = min_argmin[0];
      newLeftBound = tmp -> m_interval.geta();
    }
    tmp = tmp -> nxt;
  }

	return(newLeftBound);
}



//####### add something #######// //####### add something #######// //####### add something #######//
//####### add something #######// //####### add something #######// //####### add something #######//

void Piece::addConstant(double cst)
{
	Piece* Q_tmp = this;
	while(Q_tmp != NULL)
  {
    Q_tmp -> m_cost += cst;
    Q_tmp = Q_tmp -> nxt;
  }
}


void Piece::addPoint(Point const& pt, Robust const& robust)
{
  Piece* Q_tmp = this;


///CASE NOT ROBUST
///CASE NOT ROBUST
   if(robust.getRobustType() == "notRobust")
    { /// threshold K == INFINITY.
      while(Q_tmp != NULL)
      {
        Q_tmp -> m_cost += pt;
        Q_tmp = Q_tmp -> nxt;
      }
    }


///CASE BIWEIGHT
///CASE BIWEIGHT

  if(robust.getRobustType() == "biweight")
  { /// If a threshold K != INFINITY is present and slope a = 0
    double K = robust.getThreshold();
    ///Create an Interval new_interval with the 2 solutions of equation cost = K + minimum.
    CostGauss tempcost = CostGauss(); ///Empty cost
    tempcost += pt;   ///Adding the new Point pt
    tempcost += -K*K;

    Interval new_interval = tempcost.intervalInterRoots();

    /// Putting the bounds in variables A and B
    double AK = new_interval.geta();
    double BK = new_interval.getb();

    /// bounds of the tmp Piece
    double tmpA = 0;
    double tmpB = 0;

    int cas = 0;

    while(Q_tmp != NULL)
    {
      tmpA = Q_tmp -> m_interval.geta();
      tmpB = Q_tmp -> m_interval.getb();

      if(tmpB <= AK || BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;}///priority to AK over BK between tempA and tempB.

      switch(cas)
      {
        case 1 : Q_tmp -> m_cost += K*K;
          break;
        case 2 : Q_tmp -> m_cost += pt;
          break;
        case 3 :
        {
          ///copying tmp in new_piece
          Piece* new_piece = new Piece(Q_tmp);
          new_piece -> nxt = Q_tmp -> nxt;
              ///adding the cost to the rupt on the left
          Q_tmp -> m_cost += pt;
              ///changing interval bounds
          Q_tmp -> m_interval.setb(BK);
          new_piece -> m_interval.seta(BK);
             ///changing pointers
          Q_tmp -> nxt = new_piece;
          //delete(new_piece);
          break;
        }

        case 4 :
        {
          ///copying tmp in new_piece
          Piece* new_piece = new Piece(Q_tmp);
          new_piece -> nxt = Q_tmp -> nxt;
              ///adding the cost to the Piece on the left
          Q_tmp -> m_cost += K*K;
              ///changing interval bounds
          Q_tmp -> m_interval.setb(AK);
          new_piece -> m_interval.seta(AK);
             ///changing pointers
          Q_tmp -> nxt = new_piece;
          //delete(new_piece);
          break;
        }
      }
      Q_tmp = Q_tmp -> nxt;
    }

  }


///CASE Huber
///CASE Huber
 if(robust.getRobustType() == "Huber")
  { /// If a threshold K != INFINITY is present and slope a > 0
    double K = robust.getThreshold();

    ///Create an Interval new_interval with the 2 solutions of equation cost = K + minimum.
    CostGauss tempcost = CostGauss(); ///Empty cost
    tempcost += pt;   ///Adding the new Point pt
    tempcost += -K*K;

    Interval new_interval = tempcost.intervalInterRoots();

    /// Putting the bounds in variables A and B
    double AK = new_interval.geta();
    double BK = new_interval.getb();

    /// bounds of the tmp Piece
    double tmpA = 0;
    double tmpB = 0;

    int cas = 0;

    while(Q_tmp != NULL)
    {
      tmpA = Q_tmp -> m_interval.geta();
      tmpB = Q_tmp -> m_interval.getb();

      if(tmpB <= AK){cas = 0;}
      if(BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;}///priority to AK over BK between tempA and tempB.

      switch(cas)
      {
        case 0 : Q_tmp -> m_cost.addL1(pt, robust, -1);
          break;
        case 1 : Q_tmp -> m_cost.addL1(pt, robust, 1);
          break;
        case 2 : Q_tmp -> m_cost += pt;
          break;
        case 3 :
        {
          ///copying tmp in new_piece
          Piece* new_piece = new Piece(Q_tmp);
          new_piece -> nxt = Q_tmp -> nxt;
              ///adding the cost to the rupt on the left
          Q_tmp -> m_cost += pt;
              ///changing interval bounds
          Q_tmp -> m_interval.setb(BK);
          new_piece -> m_interval.seta(BK);
             ///changing pointers
          Q_tmp -> nxt = new_piece;
          break;
        }

        case 4 :
        {
          ///copying tmp in new_piece
          Piece* new_piece = new Piece(Q_tmp);
          new_piece -> nxt = Q_tmp -> nxt;
              ///adding the cost to the Piece on the left
          Q_tmp -> m_cost.addL1(pt, robust, -1);
              ///changing interval bounds
          Q_tmp -> m_interval.setb(AK);
          new_piece -> m_interval.seta(AK);
             ///changing pointers
          Q_tmp -> nxt = new_piece;
          break;
        }
      }
      Q_tmp = Q_tmp -> nxt;
    }
  }

}


//####### edge_constraint #######// //####### edge_constraint #######// //####### edge_constraint #######//
//####### edge_constraint #######// //####### edge_constraint #######// //####### edge_constraint #######//
//####### edge_constraint #######// //####### edge_constraint #######// //####### edge_constraint #######//
//####### edge_constraint #######// //####### edge_constraint #######// //####### edge_constraint #######//

Piece* Piece::edge_constraint(Edge const& edge, int newLabel, Bound const& bound)
{
	/// RETURN RESPONSE : first Piece address
	Piece* response;  /// first Piece to return

  /// EDGE PARAMETER
  std::string edge_ctt = edge.getConstraint();
  double edge_parameter = edge.getParameter();
  double edge_beta = edge.getBeta();
  int parentStateLabel = edge.getState1(); ///parentStateLabel = state to associate
  double m = bound.getm();
  double M = bound.getM();

  //###############################################################
  if(edge_ctt == "std" && bound.getIsConstrained() == false)
  {
    response = operator_std(newLabel, parentStateLabel, bound); ///A constant function
  }

  //###############################################################
  if(edge_ctt == "std" && bound.getIsConstrained() == true)
  {
    response = operator_stdConstr(newLabel, parentStateLabel, bound); ///A constant function
  }

  //###############################################################
  if(edge_ctt == "down")
  {
    response = operator_up(newLabel, parentStateLabel);
    if(edge_parameter > 0){response = response -> shift_left(edge_parameter, m);}
  }

  //###############################################################
  if(edge_ctt == "up")
  {
    response = operator_down(newLabel, parentStateLabel);
    if(edge_parameter > 0){response = response -> shift_right(edge_parameter, M);}
  }

  //###############################################################
  if(edge_ctt == "absInf")
  {
    response = operator_up(newLabel, parentStateLabel);
    response = response -> shift_right(edge_parameter, M);

    Piece* pieceDOWN = operator_down(newLabel, parentStateLabel);
    pieceDOWN = pieceDOWN -> shift_left(edge_parameter, m);

    response = response -> max_function(pieceDOWN, M);
    delete(pieceDOWN);
  }

  //###############################################################
  if(edge_ctt == "absSup")
  {
    response = operator_up(newLabel, parentStateLabel);
    response = response -> shift_left(edge_parameter, m);

    Piece* pieceDOWN = operator_down(newLabel, parentStateLabel);
    pieceDOWN = pieceDOWN -> shift_right(edge_parameter, M);

    response = response -> min_function(pieceDOWN, M);
    delete(pieceDOWN);
  }

  ///############  adding the constant beta (= edge_beta) ############
  response -> addConstant(edge_beta);

	return(response);
}


// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT ///
// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT /// SHIFT ///
// WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY //
// WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY // WITH NO COPY //

// shift_right // shift_right // shift_right // shift_right // shift_right // shift_right //
// shift_right // shift_right // shift_right // shift_right // shift_right // shift_right //

Piece* Piece::shift_right(double parameter, double maxi)
{
  Piece* response = this;
  Piece* Q_tmp = response;

  while(Q_tmp -> m_interval.getb() + parameter < maxi)
  {
    ///MOVE bounds
    Q_tmp -> m_interval.setb(Q_tmp -> m_interval.getb() + parameter);
    Q_tmp -> nxt -> m_interval.seta(Q_tmp -> nxt -> m_interval.geta() + parameter); ///Q_tmp -> nxt != NULL

    ///MOVE Cost
    Q_tmp -> getRefCost().shift(parameter);

    Q_tmp = Q_tmp -> nxt;
  }

  ///last Piece
  Q_tmp -> getRefCost().shift(parameter);
  Q_tmp -> m_interval.setb(maxi);

  Piece* tmp_last = Q_tmp -> nxt;
  delete(tmp_last);
  tmp_last = NULL;
  Q_tmp -> nxt = NULL;

  return(response);
}


// shift_left // shift_left // shift_left // shift_left // shift_left // shift_left //
// shift_left // shift_left // shift_left // shift_left // shift_left // shift_left //

Piece* Piece::shift_left(double parameter, double mini)
{
	Piece* Q_tmp = this;
	Piece* tmpToDelete = Q_tmp;
	Piece* ToDelete = Q_tmp;
	Piece* first_piece;

  ///Piece to delete
  while((Q_tmp -> m_interval.getb() - parameter < mini) && (Q_tmp -> nxt != NULL))
  {
    ToDelete = Q_tmp;
    Q_tmp = Q_tmp -> nxt;
  }

  ///delete the non used Pieces
  if(ToDelete != Q_tmp)
  {
    ToDelete -> nxt = NULL;
    delete(tmpToDelete);
    tmpToDelete = NULL;
  }

  ///First Piece to shift
  first_piece = Q_tmp;
  Q_tmp -> m_interval.seta(mini);

  while(Q_tmp -> nxt != NULL)
  {
    ///MOVE bounds
    Q_tmp -> m_interval.setb(Q_tmp -> m_interval.getb() - parameter);
    Q_tmp -> nxt -> m_interval.seta(Q_tmp -> nxt -> m_interval.geta() - parameter);

    ///MOVE Cost
    Q_tmp -> getRefCost().shift(-parameter);

    Q_tmp = Q_tmp -> nxt;
  }

  Q_tmp -> getRefCost().shift(-parameter);

  return(first_piece);

}



//####### operator_functions #######// //####### operator_functions #######// //####### operator_functions #######//
//####### operator_functions #######// //####### operator_functions #######// //####### operator_functions #######//
//####### operator_functions #######// //####### operator_functions #######// //####### operator_functions #######//
//####### operator_functions #######// //####### operator_functions #######// //####### operator_functions #######//




//####### operator_std_min_argmin #######// //####### operator_std_min_argmin #######// //####### operator_std_min_argmin #######//
//####### operator_std_min_argmin #######// //####### operator_std_min_argmin #######// //####### operator_std_min_argmin #######//


Piece* Piece::operator_std_min_argmin(int newLabel, int& labelmin, double& argmini, Bound const& bound)
{
  Piece* Q_tmp = this;
	Piece* Q_std = new Piece();

  /// initialization of the response piece
  Q_std -> m_interval = Interval(bound.getm(), bound.getM());
  Q_std -> m_cost = CostGauss();

  double current_min;
  double global_min = Q_tmp -> m_cost.minimum(); ///global minimum to find
  argmini = Q_tmp -> m_cost.arg_minimum(); ///argminimum to find
  labelmin = Q_tmp -> m_info.getLabel();


  Q_tmp = Q_tmp -> nxt;
  while(Q_tmp != NULL)
  {
    current_min = Q_tmp -> m_cost.minimum();  ///find min
    if(current_min < global_min)
    {
      global_min = current_min; argmini = Q_tmp -> m_cost.arg_minimum();
      labelmin = Q_tmp -> m_info.getLabel();
    }
    Q_tmp = Q_tmp -> nxt;
  }

  Q_std -> addConstant(global_min);

  /// TRACKING info
  Q_std -> m_info.setTrack(newLabel, -1, -1); ///set Track

	return(Q_std);
}



//####### operator_stdConstr_min_argmin #######// //####### operator_stdConstr_min_argmin #######// //####### operator_stdConstr_min_argmin #######//
//####### operator_stdConstr_min_argmin #######// //####### operator_stdConstr_min_argmin #######// //####### operator_stdConstr_min_argmin #######//

Piece* Piece::operator_stdConstr_min_argmin(int newLabel, int& labelmin, double& argmini, Bound const& bound)
{
  Piece* Q_tmp = this;
	Piece* Q_std = new Piece();

  /// initialization of the response piece
  Q_std -> m_interval = Interval(bound.getm(), bound.getM());
  Q_std -> m_cost = CostGauss();

  Interval inter_mM = Interval(bound.getm(), bound.getM());

  double current_min;
  double global_min = Q_tmp -> m_cost.minimum(); ///global minimum to find
  argmini = Q_tmp -> m_cost.arg_minimum(); ///argminimum to find
  if(!inter_mM.isInside(Q_tmp -> m_cost.arg_minimum()))
  {
    if(Q_tmp -> m_cost.arg_minimum() < inter_mM.geta()){global_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.geta()); argmini = Q_tmp -> m_interval.geta();}
    if(Q_tmp -> m_cost.arg_minimum() > inter_mM.getb()){global_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.getb()); argmini = Q_tmp -> m_interval.getb();}
  }
  labelmin = Q_tmp -> m_info.getLabel();


  Q_tmp = Q_tmp -> nxt;
  while(Q_tmp != NULL)
  {
    current_min = Q_tmp -> m_cost.minimum();  ///find min
    if(!inter_mM.isInside(Q_tmp -> m_cost.arg_minimum()))
    {
      if(Q_tmp -> m_cost.arg_minimum() < inter_mM.geta()){current_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.geta());}
      if(Q_tmp -> m_cost.arg_minimum() > inter_mM.getb()){current_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.getb());}
    }
    if(current_min < global_min)
    {
      global_min = current_min;
      argmini = Q_tmp -> m_cost.arg_minimum();
      if(Q_tmp -> m_cost.arg_minimum() < inter_mM.geta()){argmini = Q_tmp -> m_interval.geta();}
      if(Q_tmp -> m_cost.arg_minimum() > inter_mM.getb()){argmini = Q_tmp -> m_interval.getb();}
      labelmin = Q_tmp -> m_info.getLabel();
    }
    Q_tmp = Q_tmp -> nxt;
  }

  Q_std -> addConstant(global_min);

  /// TRACKING info
  Q_std -> m_info.setTrack(newLabel, -1, -1); ///set Track

	return(Q_std);

}


//####### operator std #######// //####### operator std #######// //####### operator std #######//
//####### operator std #######// //####### operator std #######// //####### operator std #######//


Piece* Piece::operator_std(int newLabel, int parentStateLabel, Bound const& bound)
{
  Piece* Q_tmp = this;
	Piece* Q_std = new Piece();

  /// initialization of the response piece
  Q_std -> m_interval = Interval(bound.getm(), bound.getM());
  Q_std -> m_cost = CostGauss();

  double global_min = Q_tmp -> m_cost.minimum(); ///global minimum to find
  double current_min;
  int parentPosition = 1;
  int counter = 1;
  Q_tmp = Q_tmp -> nxt; counter = counter + 1;

  while(Q_tmp != NULL)
  {
    current_min = Q_tmp -> m_cost.minimum();  ///find min
    if(current_min < global_min){global_min = current_min; parentPosition = counter;}
    Q_tmp = Q_tmp -> nxt; counter = counter + 1;
  }

  Q_std -> addConstant(global_min);

  /// TRACKING info
  Q_std -> m_info.setTrack(newLabel, parentStateLabel, parentPosition); ///set Track

	return(Q_std);
}

//####### operator stdConstr #######// //####### operator stdConstr #######// //####### operator stdConstr #######//
//####### operator stdConstr #######// //####### operator stdConstr #######// //####### operator stdConstr #######//


Piece* Piece::operator_stdConstr(int newLabel, int parentStateLabel, Bound const& bound)
{
  Piece* Q_tmp = this;
	Piece* Q_std = new Piece();

  /// initialization of the response piece
  Q_std -> m_interval = Interval(bound.getm(), bound.getM());
  Q_std -> m_cost = CostGauss();

  double global_min; ///global minimum to find
  double current_min;
  int parentPosition = 1;
  int counter = 1;

  Interval inter_mM = Interval(bound.getm(), bound.getM());

  global_min = Q_tmp -> m_cost.minimum(); ///global minimum to find
  if(!inter_mM.isInside(Q_tmp -> m_cost.arg_minimum()))
    {
      if(Q_tmp -> m_cost.arg_minimum() < inter_mM.geta()){global_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.geta());}
      if(Q_tmp -> m_cost.arg_minimum() > inter_mM.getb()){global_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.getb());}
    }
  Q_tmp = Q_tmp -> nxt; counter = counter + 1;

  while(Q_tmp != NULL)
  {
    current_min = Q_tmp -> m_cost.minimum();  ///find min
    if(!inter_mM.isInside(Q_tmp -> m_cost.arg_minimum()))
    {
      if(Q_tmp -> m_cost.arg_minimum() < inter_mM.geta()){current_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.geta());}
      if(Q_tmp -> m_cost.arg_minimum() > inter_mM.getb()){current_min = Q_tmp -> m_cost.point_eval(Q_tmp -> m_interval.getb());}
    }

    if(current_min < global_min){global_min = current_min; parentPosition = counter;}
    Q_tmp = Q_tmp -> nxt; counter = counter + 1;
  }

  Q_std -> addConstant(global_min);

  /// TRACKING info
  Q_std -> m_info.setTrack(newLabel, parentStateLabel, parentPosition); ///set Track

	return(Q_std);
}



//####### operator_up #######// //####### operator_up #######// //####### operator_up #######//
//####### operator_up #######// //####### operator_up #######// //####### operator_up #######//

Piece* Piece::operator_up(int newLabel, int parentStateLabel)
{
  int length = 0;
  Piece* Q_up;
	Piece* Q_up_copy = this -> copy(length);
	//std::cout << length << std::endl;
	Q_up_copy = Q_up_copy -> reverse();
  Q_up = Q_up_copy -> operator_down(newLabel, parentStateLabel);
  delete(Q_up_copy);
  Q_up = Q_up -> reverse(length);
	return(Q_up);
}


//####### operator_down #######// //####### operator_down #######// //####### operator_down #######//
//####### operator_down #######// //####### operator_down #######// //####### operator_down #######//

Piece* Piece::operator_down(int newLabel, int parentStateLabel)
{
  Piece* Q_tmp = this;
  Piece* BUILD = new Piece();

    //std::cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV" << std::endl;
    //Q_tmp -> show();

  int counter = 1;

  Track trackDown = Track(newLabel, parentStateLabel, counter);
  BUILD -> m_info.setTrack(trackDown); ///set Track

  Piece* Q_down = BUILD; ///response Piece

  //std::cout << BUILD << std::endl;
  ///First Piece
  double bound = Q_tmp -> m_interval.geta();
  BUILD -> m_interval.seta(bound);
  BUILD -> m_interval.setb(bound);
  //BUILD -> m_cost = CostGauss(); => in new Piece();
  double currentPoint = Q_tmp -> m_cost.point_eval(bound);
  BUILD -> m_cost += currentPoint;

  ///First Piece bool constPiece = is the first Piece constant? If cost increasing at bound, constPiece = true
  bool constPiece = false;
  if(Q_tmp -> m_cost.arg_minimum() < bound){constPiece = true;}


  Interval decreasingInterval = Interval(); /// for interval building

  while(Q_tmp != NULL)
  {
    ///Interval on Q_tmp to create
    decreasingInterval = Q_tmp -> m_cost.intervalMinLess(currentPoint, bound, constPiece); ///"decreasing" interval
    decreasingInterval = decreasingInterval.intersection(Q_tmp -> getInterval()); ///decreasingInterval = intersection of decreasingInterval (=intervalMinLess) and interval of  Q_tmp

    //////////////////////////////////////////////////////////////////////////////////////////
    if(decreasingInterval.isEmpty() == false){trackDown.setPosition(counter);}
    BUILD = BUILD -> pastePiece(Q_tmp, decreasingInterval, trackDown); ///add new Piece to BUILD

    //std::cout << BUILD << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////

    bound = BUILD -> m_interval.getb(); ///new rightBound
    currentPoint = BUILD -> m_cost.point_eval(bound); ///new currentPoint (=the minimum)

    ///update constPiece (true => false => true is possible)
    if(constPiece == true){if(decreasingInterval.isEmpty() == false){constPiece = false;}}
    if(constPiece == false)
    {
      if(decreasingInterval.getb() < Q_tmp -> getInterval().getb()){constPiece = true;}
    }

    Q_tmp = Q_tmp -> nxt;
    counter = counter + 1;
  }

  //std::cout << BUILD << std::endl;
  //std::cout << BUILD -> nxt << std::endl;
  BUILD = BUILD -> nxt; // =NULL
  delete(BUILD);

  return(Q_down);
}



//####### pasteConstraint #######// //####### pasteConstraint #######// //####### pasteConstraint #######//
//####### pasteConstraint #######// //####### pasteConstraint #######// //####### pasteConstraint #######//


Piece* Piece::pastePiece(const Piece* Q, Interval const& decrInter, Track const& newTrack)
{
  //std::cout << "pastePiece" << std::endl;
  //std::cout << this << std::endl;
  Piece* BUILD = this;
  //std::cout << BUILD << std::endl;

  /// decreasingInterval = (a,b)
  /// Q -> m_interval = (m_a,m_b)

  if(decrInter.isEmpty())
  {
    BUILD -> m_interval.setb(Q -> getInterval().getb());
    //std::cout << "1) D==0 : " << BUILD -> getTrack().getPosition() << std::endl;
  }
  else
  {
    BUILD -> m_interval.setb(decrInter.geta());    ///if a > m_a, no change otherwise

    ///ADD of the PIECE Q (troncated)
    if(BUILD -> m_interval.isEmpty()) ///if BUILD empty
    {
      BUILD -> m_interval.setb(decrInter.getb());
      BUILD -> m_cost = Q -> m_cost;
      BUILD -> m_info.setTrack(newTrack);
      //std::cout << "2) D!=0 & enterI == 0 : " << BUILD -> getTrack().getPosition() << std::endl;

    }
    else
    {
      Piece* NewQ = new Piece(newTrack, decrInter, Q -> m_cost);
      BUILD -> nxt = NewQ;
      BUILD = NewQ;
      //++++POSITION OF Q
      //std::cout << "3) D!=0 & enterI != 0 : " << BUILD -> getTrack().getPosition() << std::endl;
    }

    if(!((Q -> nxt == NULL) && (decrInter.getb() == Q -> m_interval.getb())))
    {
      ///OUTPUT PIECE. value = const, Interval = (decrInter.getb(), Q -> m_interval.getb()) = or != empty
      double outputValue = Q -> m_cost.point_eval(decrInter.getb());
      Piece* PieceOut = new Piece(newTrack, Interval(decrInter.getb(), Q -> m_interval.getb()), CostGauss());
      PieceOut -> m_cost += outputValue;
      BUILD -> nxt = PieceOut;
      BUILD = PieceOut;
      //std::cout << "     D!=0 & Db == Qb : " << BUILD -> getTrack().getPosition() << " & " << decrInter.getb() << " - " << Q -> m_interval.getb() << std::endl;
    }


  }

  //std::cout << BUILD << std::endl;

  return(BUILD);

}




// COMPARISON OF FUNCIONS /// COMPARISON OF INTERNAL FUNCIONS /// COMPARISON OF INTERNAL FUNCIONS
// COMPARISON OF FUNCIONS /// COMPARISON OF INTERNAL FUNCIONS /// COMPARISON OF INTERNAL FUNCIONS
// min_function /// min_function /// min_function /// min_function /// min_function ///
// min_function /// min_function /// min_function /// min_function /// min_function ///

Piece* Piece::min_function(Piece* Q2, double M)
{
  //"MOST OF THE TIME" : Q2 > Q1
	Piece* Q1 = this;  /// first Piece to compare
	Piece* Q12 = new Piece();
  Piece* response = Q12;

  ///
  /// SET Q2_Minus_Q1
  ///
  double firstLeftBound = Q1 -> m_interval.geta();

  /// Q2_Minus_Q1 -> = 1 if Q2 >= Q1 at the current point, - 1 if Q2 < Q1
  ///
  /// Q12 first Piece
  ///
  Q12 -> m_interval = Interval(firstLeftBound,firstLeftBound);

  int Bound_Q2_Minus_Q1 = 0;

  while(Q1 != NULL)
  {
    ///Q12 = Piece with an interval but no cost no label
    /// Bound_Q2_Minus_Q1
    /// = 0 if bound interval Q2 - bound interval Q1 == 0 : Q1 and Q2 stop
    /// = 1 if bound interval Q2 - bound interval Q1 > 0: Q1 stops
    /// = -1 if bound interval Q2 - bound interval Q1 < 0 : Q2 stops
    Bound_Q2_Minus_Q1 = -1;

    while(Bound_Q2_Minus_Q1 == -1)
    {
      /// right bound
      if(Q2 -> m_interval.getb() - Q1 -> m_interval.getb() > 0){Bound_Q2_Minus_Q1 = 1;}
      if(Q1 -> m_interval.getb() - Q2 -> m_interval.getb() == 0){Bound_Q2_Minus_Q1 = 0;}

      Q12 = Q12 -> pieceGenerator(Q1, Q2, Bound_Q2_Minus_Q1, M); ///add new Piece(s) to Q12

      if(Bound_Q2_Minus_Q1 < 1){Q2 = Q2 -> nxt;}
    }

  Q1 = Q1 -> nxt;
  }

  delete(this);
	return(response);
}


// max_function /// max_function /// max_function /// max_function /// max_function ///
// max_function /// max_function /// max_function /// max_function /// max_function ///

Piece* Piece::max_function(Piece* Q2, double M)
{
	Piece* Q1 = this;  /// first Piece to compare

	Q1 -> opposition();
	Q2 -> opposition();

  Piece* Q12 = Q1 -> min_function(Q2, M);
	Q12 -> opposition();

	///Q1 -> opposition(); Q1 already deleted
	Q2 -> opposition();

	return(Q12);
}


//####### pieceGenerator #######// //####### pieceGenerator #######// //####### pieceGenerator #######//
//####### pieceGenerator #######// //####### pieceGenerator #######// //####### pieceGenerator #######//

Piece* Piece::pieceGenerator(Piece* Q1, Piece* Q2, int Bound_Q2_Minus_Q1, double M)
{
  Piece* BUILD = this;

  //// INFORMATION interToPaste
  // Interval interToPaste = interval on which we build BUILD
  // interToPaste : right BUILD -> min(right Q1, right Q2)
  Interval interToPaste;
  interToPaste.seta(BUILD -> m_interval.getb()); ///enter Piece BUILD :
  if(Bound_Q2_Minus_Q1 == -1){interToPaste.setb(Q2 -> m_interval.getb());}else{interToPaste.setb(Q1 -> m_interval.getb());}

  //// INFORMATION interRoots
  // Interval interRoots (Q1 - Q2)
  Interval interRoots = (Q1 -> m_cost.minus(Q2 -> m_cost)).intervalInterRoots();

  //// INFORMATION change
  // int change = 0, 1 or 2 change-points
  // if bounds of interRoots close to bounds of interToPaste > 1e-12
  int change = 0;
  if((interRoots.geta() > interToPaste.geta() + 1e-12)&&(interRoots.geta() + 1e-12 < interToPaste.getb())){change = change + 1;}
  if((interRoots.getb() > interToPaste.geta() + 1e-12)&&(interRoots.getb() + 1e-12 < interToPaste.getb())){change = change + 1;}

  ///Security steps: length interRoots very small < 1e-4 DANGER DANGER DANGER DANGER if put to < 1e-12 ???
  if(interRoots.getb() - interRoots.geta() < 1e-4)
  {
    change = 0;
    interRoots.seta(INFINITY);
    interRoots.setb(INFINITY);
  }
  //std::cout << "change: " << change << std::endl;

  //// CONSTRUCTION
  // CONSTRUCTION
  double centerPoint;
  int Q2_Minus_Q1;  ///Sign of Q2 - Q1
  CostGauss testCost;
  bool test;

  //std::cout << "change: " << change << std::endl;

  switch(change)
  {
    /// IF WE ADD 0 PIECE
    case 0 :
    {
      /// Possible inversion => test Q2_Minus_Q1 at centerPoint
      centerPoint = interToPaste.internPoint();
      Q2_Minus_Q1 = (Q1 -> m_cost).sign_Q2_Minus_Q1(Q2 -> m_cost, centerPoint);

      if (BUILD -> getInterval().isEmpty() == true) /// IF BUILD interval = empty
      {
        //PROLONGATION
        BUILD -> m_interval.setb(interToPaste.getb());
        if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
        if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}
      }
      else
      {
        ///SECURITY step
        testCost = BUILD -> getCost();
        if(Q2_Minus_Q1 == 1){test = Q1 -> getCost().isEqual(testCost);}
        if(Q2_Minus_Q1 == -1){test = Q2 -> getCost().isEqual(testCost);}

        if (test == true) ///no pb with the cost
        {
          //PROLONGATION
          BUILD -> m_interval.setb(interToPaste.getb());
          if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
          if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}
        }
        else ///pb with the cost -> we stop BUILD interval at interToPaste left -> we create a new piece
        {
          //CONSTRUCTION newPiece
          BUILD -> m_interval.setb(interToPaste.geta());
          Piece* newPiece = new Piece();
          newPiece -> m_interval = interToPaste;
          if(Q2_Minus_Q1 == 1){newPiece -> m_cost = Q1 -> m_cost; newPiece -> m_info = Q1 -> m_info;}
          if(Q2_Minus_Q1 == -1){newPiece -> m_cost = Q2 -> m_cost; newPiece -> m_info = Q2 -> m_info;}
          BUILD -> nxt = newPiece;
          BUILD = newPiece;
        }
      }
    break;
    }

    /// IF WE ADD 1 PIECE
    case 1 :
    {
      //PROLONGATION theChangePoint
      double theChangePoint;
      if(interToPaste.geta() < interRoots.geta()){theChangePoint = interRoots.geta();}else{theChangePoint = interRoots.getb();}

      //// FIND the winner on the new piece
      // centerPoint = centre (left interToPaste, right theChangePoint)
      centerPoint = Interval(interToPaste.geta(), theChangePoint).internPoint();
      Q2_Minus_Q1 = (Q1 -> m_cost).sign_Q2_Minus_Q1(Q2 -> m_cost, centerPoint);

      if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}

      BUILD -> m_interval.setb(theChangePoint);

      //std::cout << centerPoint <<  " &&& " << theChangePoint<< std::endl;
      //std::cout << "interRoots "; interRoots.show();
      //std::cout << "interToPaste "; interToPaste.show();
      //std::cout << "BUILD "; BUILD -> m_interval.show();
      //std::cout << "Q1 "; Q1 -> m_interval.show();
      //std::cout << "Q2 "; Q2 -> m_interval.show();

      //CONSTRUCTION newPiece
      Piece* newPiece = new Piece();
      newPiece -> m_interval = Interval(theChangePoint, interToPaste.getb());

      centerPoint = newPiece -> getInterval().internPoint();
      Q2_Minus_Q1 = (Q1 -> m_cost).sign_Q2_Minus_Q1(Q2 -> m_cost, centerPoint);

      if(Q2_Minus_Q1 == 1){newPiece -> m_cost = Q1 -> m_cost; newPiece -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){newPiece -> m_cost = Q2 -> m_cost; newPiece -> m_info = Q2 -> m_info;}
      BUILD -> nxt = newPiece;
      BUILD = newPiece;

      break;
    }

    /// IF WE ADD 2 PIECES
    case 2 :
    {
      //PROLONGATION theChangePoint
      //FIND the winner on the newpiece1 (the central piece defined on interval interRoots)
      centerPoint = interRoots.internPoint();
      Q2_Minus_Q1 = -(Q1 -> m_cost).sign_Q2_Minus_Q1(Q2 -> m_cost, centerPoint); ///INVERSION!!!

      if(Q2_Minus_Q1 == 1){BUILD -> m_cost = Q1 -> m_cost; BUILD -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){BUILD -> m_cost = Q2 -> m_cost; BUILD -> m_info = Q2 -> m_info;}

      BUILD -> m_interval.setb(interRoots.geta());

      //CONSTRUCTION newPiece1
      Q2_Minus_Q1 = -Q2_Minus_Q1; ///INVERSION!!!
      Piece* newPiece1 = new Piece();
      newPiece1 -> m_interval = interRoots;
      if(Q2_Minus_Q1 == 1){newPiece1 -> m_cost = Q1 -> m_cost; newPiece1 -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){newPiece1 -> m_cost = Q2 -> m_cost; newPiece1 -> m_info = Q2 -> m_info;}
      BUILD -> nxt = newPiece1;
      BUILD = newPiece1;

      Q2_Minus_Q1 = -Q2_Minus_Q1;
      //CONSTRUCTION newPiece2
      Piece* newPiece2 = new Piece();
      newPiece2 -> m_interval = Interval(interRoots.getb(), interToPaste.getb());
      if(Q2_Minus_Q1 == 1){newPiece2 -> m_cost = Q1 -> m_cost; newPiece2 -> m_info = Q1 -> m_info;}
      if(Q2_Minus_Q1 == -1){newPiece2 -> m_cost = Q2 -> m_cost; newPiece2 -> m_info = Q2 -> m_info;}
      BUILD -> nxt = newPiece2;
      BUILD = newPiece2;
      break;
    }
  }

  //CONSTRUCTION outPiece
  ///Need of a last "OUT" Piece if we have reached the end of the Piece Q1 or Piece Q2
  if((((Q2_Minus_Q1 == 1)&&(Bound_Q2_Minus_Q1 >= 0)) || ((Q2_Minus_Q1 == -1)&&(Bound_Q2_Minus_Q1 <= 0))) && (interToPaste.getb() != M))
  {
    Piece* outPiece = new Piece();
    outPiece -> m_interval = Interval(interToPaste.getb(), interToPaste.getb());
    BUILD -> nxt = outPiece;
    BUILD = outPiece;
  }

  return(BUILD);
}



//####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######//
//####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######//
///We test all the Piece

std::vector<double> Piece::get_min_argmin_label_state_position_final()
{
  std::vector<double> response(5,0);
  Piece* Q_tmp = this;

///INITIALIZATION
  std::vector<double> current_min_argmin = Q_tmp -> Piece_min_argmin();
  double global_min = current_min_argmin[0]; ///global minimum to find
  double global_argmin = current_min_argmin[1]; ///global argminimum to find

  int label = Q_tmp -> m_info.getLabel();
  int state = Q_tmp -> m_info.getState();
  int position = Q_tmp -> m_info.getPosition();

  Q_tmp = Q_tmp -> nxt;

///LOOP TESTS
  while(Q_tmp != NULL)
  {
    current_min_argmin = Q_tmp -> Piece_min_argmin();
    if(current_min_argmin[0] < global_min)
    {
      global_min = current_min_argmin[0];
      global_argmin = current_min_argmin[1];
      label = Q_tmp -> m_info.getLabel();
      state = Q_tmp -> m_info.getState();
      position = Q_tmp -> m_info.getPosition();
    }
    Q_tmp = Q_tmp -> nxt;
  }

  response[0] = global_min;
  response[1] = global_argmin;
  response[2] = label;
  response[3] = state;
  response[4] = position;

  return(response);
}



//####### min_argmin_label_state_position #######// //####### min_argmin_label_state_position #######// //####### min_argmin_label_state_position #######//
//####### min_argmin_label_state_position #######// //####### min_argmin_label_state_position #######// //####### min_argmin_label_state_position #######//
///We test only the piece at position i


std::vector<double> Piece::get_min_argmin_label_state_position(int i, Interval const& constrainedInter, bool out, bool& forced, bool isBoundConstrained)
{
  std::vector<double> response(5,0);
  std::vector<double> min_argmin;
  Piece* Q_tmp = this;

  ///Go to the right Piece
  int pos = 1;
  while(pos != i){Q_tmp = Q_tmp -> nxt; pos = pos + 1;}

  //std::cout << "ZZZZZZZZZZZZZZZZZZ" << std::endl;
  //Q_tmp->showOne();
  //std::cout << "ZZZZZZZZZZZZZZZZZZ" << std::endl;

  min_argmin = Q_tmp -> Piece_min_argmin();
  response[0] = min_argmin[0];
  response[1] = min_argmin[1];
  response[2] = Q_tmp -> m_info.getLabel();
  response[3] = Q_tmp -> m_info.getState();
  response[4] = Q_tmp -> m_info.getPosition();

  ///force arg_min (mean in segmentation) to fit the constaints
  /// out = false : we want an argmin in the interval constrainedInter
  /// out = true : we want an argmin outside the interval constrainedInter
  if(out == false)
  {
    if(constrainedInter.geta() >= response[1]){response[1] = constrainedInter.geta(); forced = true;}
    if(constrainedInter.getb() <= response[1]){response[1] = constrainedInter.getb(); forced = true;}
  }

  if(out == true)
  {
    if((constrainedInter.geta() < response[1]) && (response[1] < constrainedInter.getb()))
    {
    forced = true;
    if(response[1] - constrainedInter.geta() < constrainedInter.getb() - response[1]){response[1] = constrainedInter.geta();}
      else{response[1] = constrainedInter.getb();}
    }
  }
  return(response);
}



//####### min_argmin_label #######// //####### min_argmin_label #######// //####### min_argmin_label #######//
//####### min_argmin_label #######// //####### min_argmin_label #######// //####### min_argmin_label #######//

std::vector<double> Piece::get_min_argmin_label(double rightBound, bool& forced, bool isBoundConstrained)
{
  std::vector<double> response(3,0);
  Piece* Q_tmp = this;

///INITIALIZATION
  std::vector<double> current_min_argmin = Q_tmp -> Piece_min_argmin();
  double global_min = current_min_argmin[0]; ///global minimum to find
  double global_argmin = current_min_argmin[1]; ///global argminimum to find
  int label = Q_tmp -> m_info.getLabel();


  if(rightBound <= Q_tmp -> getInterval().geta())
  {
    global_min = Q_tmp -> getCost().point_eval(Q_tmp -> getInterval().geta());
    global_argmin = Q_tmp -> getInterval().geta();
    forced = true;
  }
  if(global_argmin >= rightBound)
  {
    global_min = Q_tmp -> getCost().point_eval(rightBound);
    global_argmin = rightBound;
    forced = true;
  }
  //std::cout<< global_argmin << " --- " << label << std::endl;


  Q_tmp = Q_tmp -> nxt;

///LOOP TESTS
  while((Q_tmp != NULL) && (Q_tmp -> getInterval().geta() <= rightBound))
  {
  current_min_argmin = Q_tmp -> Piece_min_argmin();
    if(current_min_argmin[0] < global_min)
    {
      global_min = current_min_argmin[0];
      global_argmin = current_min_argmin[1];
      label = Q_tmp -> m_info.getLabel();

      if(global_argmin >= rightBound)
      {
        //std::cout<<"azazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazazaz"<<std::endl;
        //std::cout<< rightBound <<std::endl;
        //std::cout<< "new" << global_argmin <<std::endl;

        global_min = Q_tmp -> getCost().point_eval(rightBound);
        global_argmin = rightBound;
        forced = true;
      }
    }
    Q_tmp = Q_tmp -> nxt;
  }

  response[0] = global_min;
  response[1] = global_argmin;
  response[2] = label;

  return(response);
}






// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE ///
// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE /// SAVE ///

void Piece::save(std::ostream &flux)
{
	Piece* tmp = this;
	flux << "L" << " " <<  "a" << " " << "b" << " " << "A" << " " << "B" << " " << "C" << std::endl;
	while(tmp != NULL)
  {
    flux << tmp -> m_info.getLabel() << " "<< tmp -> m_interval.geta() << " " << tmp -> m_interval.getb() << " " << tmp -> m_cost.getM_A() <<  " " << tmp -> m_cost.getM_B() << " " << tmp -> m_cost.getConstant() << std::endl;
    tmp = tmp -> nxt;
  }
}



void Piece::show()
{
  int i = 1;
	Piece* tmp = this;
	//if(tmp == NULL){std::cout << termcolor::red << "#NULL EMPTY POINTER# "<< termcolor::reset << std::endl;}
	while(tmp != NULL)
  {
	  //std::cout << i << "#";
	  //std::cout << tmp;
	  //std::cout << termcolor::cyan << "#LABEL# "<< tmp -> m_info.getLabel() << " #STATE# " <<  tmp -> m_info.getState() << " POSITION " <<  tmp -> m_info.getPosition() << " " << termcolor::reset;
	  //std::cout << "#INTERVAL# "<< tmp -> m_interval.geta() << " -- " << tmp -> m_interval.getb()<< " ";
	  //tmp -> m_cost.show(tmp -> m_interval);
    tmp = tmp -> nxt;
    i = i + 1;
  }
}


void Piece::showOne()
{
	Piece* tmp = this;
  //if(tmp == NULL){std::cout << termcolor::red << "#NULL EMPTY POINTER# "<< termcolor::reset << std::endl;}
  //else
  //{
  //  std::cout << tmp;
  //  std::cout << termcolor::cyan << "#LABEL# "<< tmp -> m_info.getLabel() << " #STATE# " <<  tmp -> m_info.getState() << " POSITION " <<  tmp -> m_info.getPosition() << " " << termcolor::reset;
  //  std::cout << "#INTERVAL# "<< tmp -> m_interval.geta() << " -- " << tmp -> m_interval.getb()<< " ";
  //  tmp -> m_cost.show(tmp -> m_interval);
  //}
}


std::ostream &operator>>(std::ostream &flux, Piece* piece)
{
  piece -> save(flux);
  return(flux);
}


