/* -*- compile-command: "R CMD INSTALL .." -*- */
#include "ListPiece.h"

#include "Piece.h"
#include <iostream>
#include "stdlib.h"


ListPiece::ListPiece()
{
  head = NULL;
  currentPiece = NULL;
  lastPiece = NULL;
}

ListPiece::~ListPiece()
{
  delete(head);
  head = NULL;
  currentPiece = NULL;
  lastPiece = NULL;
}

//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///
//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///

void ListPiece::setUniquePieceCostToInfinity()
{
  head -> m_cost.constant = INFINITY;
}


//##### setNewBounds #####//////##### setNewBounds #####//////##### setNewBounds #####///
//##### setNewBounds #####//////##### setNewBounds #####//////##### setNewBounds #####///

void ListPiece::setNewBounds(Interval newBounds)
{
  double a = newBounds.geta();
  double b = newBounds.getb();
  Piece* tmp;

  //left bound
  if(a <= head -> m_interval.geta()){head -> m_interval.seta(a);}
  else
  {
    while(a > head -> m_interval.getb())
    {
      tmp = head;
      head = head -> nxt;
      tmp -> nxt = NULL;
      delete(tmp);
    }
    head -> m_interval.seta(a);
  }

  //right bound
  if(b >= lastPiece -> m_interval.getb()){lastPiece -> m_interval.setb(b);}
  else
  {
    tmp = head;
    while(tmp -> m_interval.getb() < b){tmp = tmp -> nxt;}
    tmp -> m_interval.setb(b);
    if(tmp -> nxt != NULL){delete(tmp -> nxt); tmp -> nxt = NULL;}
    lastPiece = tmp;
  }
}

//##### reset #####//////##### reset #####//////##### reset #####///
//##### reset #####//////##### reset #####//////##### reset #####///

void ListPiece::reset()
{
  while(head != NULL)
  {
    Piece* pieceToDelete = head;
    head = head -> nxt;
    pieceToDelete -> nxt = NULL;
    delete(pieceToDelete);
  }
  currentPiece = NULL;
  lastPiece = NULL;
}

//##### copy #####//////##### copy #####//////##### copy #####///
//##### copy #####//////##### copy #####//////##### copy #####///

void ListPiece::copy(ListPiece const& LP_edge)
{
  Piece* tmp = LP_edge.head;
  head = tmp -> copy();
  currentPiece = head;
  tmp = tmp -> nxt;

  while(tmp != NULL)
  {
    currentPiece -> nxt = tmp -> copy(); //copy content in Piece
    currentPiece = currentPiece -> nxt; //copy nxt pointer
    tmp = tmp -> nxt;
  }
  lastPiece = currentPiece;
}


//##### reverseAndCount #####//////##### reverseAndCount #####//////##### reverseAndCount #####///
//##### reverseAndCount #####//////##### reverseAndCount #####//////##### reverseAndCount #####///

void ListPiece::reverseAndCount(unsigned int& length)
{
  lastPiece = head;

  Piece* prev =  NULL;
  Piece* current = head;
  Piece* next = current;

  while(current != NULL)
  {
    next  = current -> nxt;
    current -> nxt = prev; /// new nxt
    prev = current;
    length = length + 1; ///add +1 to length
    current = next;
  }

  head = prev;
  currentPiece = head;
}


//##### reverseAndSetTrackPosition #####//////##### reverseAndSetTrackPosition #####//////##### reverseAndSetTrackPosition #####///
//##### reverseAndSetTrackPosition #####//////##### reverseAndSetTrackPosition #####//////##### reverseAndSetTrackPosition #####///

void ListPiece::reverseAndSetTrackPosition(unsigned int length)
{
  lastPiece = head;

  Piece* prev =  NULL;
  Piece* tmp = head;
  Piece* next = tmp;

  while(tmp != NULL)
  {
    next  = tmp -> nxt;
    tmp -> nxt = prev; /// new nxt
    prev = tmp;
    tmp -> m_info.reversePosition(length);
    tmp = next;
  }

  head = prev;
  currentPiece = head;
}


//##### addCurrentPiecePlus1NotMove #####//////##### addCurrentPiecePlus1NotMove #####//////##### addCurrentPiecePlus1NotMove #####///
//##### addCurrentPiecePlus1NotMove #####//////##### addCurrentPiecePlus1NotMove #####//////##### addCurrentPiecePlus1NotMove #####///

void ListPiece::addCurrentPiecePlus1NotMove(Piece* newPiece)
{
  newPiece -> nxt = currentPiece -> nxt;
  currentPiece -> nxt = newPiece;
}

//##### addFirstPiece #####//////##### addFirstPiece #####//////##### addFirstPiece #####///
//##### addFirstPiece #####//////##### addFirstPiece #####//////##### addFirstPiece #####///

void ListPiece::addFirstPiece(Piece* newPiece)
{
  head = newPiece;
  currentPiece = newPiece;
  lastPiece = newPiece;
}

//##### move #####//////##### move #####//////##### move #####///
//##### move #####//////##### move #####//////##### move #####///

void ListPiece::move()
{
  currentPiece = currentPiece -> nxt;
}

//##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####///
//##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####//////##### initializeCurrentPiece #####///

void ListPiece::initializeCurrentPiece()
{
  currentPiece = head;
}

//##### initializeHeadWithFirstPoint #####//////##### initializeHeadWithFirstPoint #####//////##### initializeHeadWithFirstPoint #####///
//##### initializeHeadWithFirstPoint #####//////##### initializeHeadWithFirstPoint #####//////##### initializeHeadWithFirstPoint #####///

void ListPiece::initializeHeadWithFirstPoint(Point const& pt)
{
  double* coeff = cost_coeff(pt);
  Cost costPt = Cost(coeff);
  head -> addCostAndPenalty(costPt, 0);
  delete[] coeff;
}

//##### shift #####//////##### shift #####//////##### shift #####///
//##### shift #####//////##### shift #####//////##### shift #####///

void ListPiece::shift(double parameter)
{
  Interval inter;
  initializeCurrentPiece();
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interShift(inter.geta(), parameter));
    currentPiece -> m_interval.setb(cost_interShift(inter.getb(), parameter));
    ///MOVE Cost
    cost_shift(currentPiece -> m_cost, parameter);
    move();
  }
}


//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///
//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///

void ListPiece::expDecay(double gamma)
{
  Interval inter;
  initializeCurrentPiece();
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interExpDecay(inter.geta(), gamma));
    currentPiece -> m_interval.setb(cost_interExpDecay(inter.getb(), gamma));
    ///MOVE Cost
    cost_expDecay(currentPiece -> m_cost, gamma);
    move();
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///
//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///

void ListPiece::LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel)
{
  reset(); /// build a new LP_edges from scratch

  /// only 4 types of edges : null, std, up, down
  /// EDGE PARAMETERS
  std::string edge_ctt = edge.getConstraint();
  double edge_parameter = edge.getParameter(); /// always positive
  unsigned int parentState = edge.getState1(); ///parentState = state to associate

  //################
  if(edge_ctt == "null") /// Simple copy of LP_state
  {
    copy(LP_state);
    if(edge_parameter < 1){expDecay(edge_parameter);} ///edge_parameter = exponential decay
  }

  //################
  if(edge_ctt == "std")
  {
    ///variable definition
    Piece* tmp = LP_state.head;
    double globalMin = INFINITY;
    unsigned int positionMin = 0;
    unsigned int currentCounter = 0;
    double currentMin;


    ///find the minimum
    while(tmp != NULL)
    {
      currentCounter = currentCounter + 1;
      currentMin = cost_minInterval(tmp -> m_cost, tmp -> m_interval);
      if(currentMin < globalMin){globalMin = currentMin; positionMin = currentCounter;}
      tmp = tmp -> nxt;
    }

    ///add onePiece to LP_edges
    Piece* onePiece = new Piece();
    onePiece -> m_info = Track(newLabel, parentState, positionMin);
    onePiece -> m_interval = Interval(LP_state.head -> m_interval.geta(), LP_state.lastPiece -> m_interval.getb());
    onePiece -> addCostAndPenalty(Cost(), globalMin); /// Cost() = 0
    addFirstPiece(onePiece);
  }


  //################
  if(edge_ctt == "up")
  {
    operatorUp(LP_state, newLabel, parentState);
    if(edge_parameter > 0){shift(edge_parameter);} ///edge_parameter = right decay
  }

  //################
  if(edge_ctt == "down")
  {
    ListPiece LP_stateCopy = ListPiece(); ///LP_state working copy
    LP_stateCopy.copy(LP_state);
    unsigned int length = 0;
    LP_stateCopy.reverseAndCount(length); ///reverse LP_stateCopy

    operatorDw(LP_stateCopy, newLabel, parentState); ///down operations
    reverseAndSetTrackPosition(length); ///reverse result

    if(edge_parameter > 0){shift(-edge_parameter);} ///edge_parameter = left decay
  }

}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void ListPiece::LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt)
{
  /// get edge data ///
  double K = edge.getKK();
  double a = edge.getAA();
  double edge_beta = edge.getBeta();
  /// get pt cost ///
  double* coeff = cost_coeff(pt);
  Cost costPt = Cost(coeff);

  initializeCurrentPiece();

  ///////////////////// CASE K == INF /////////////////////
  if(K == INFINITY)
  {
    while(currentPiece != NULL)
    {
      currentPiece -> addCostAndPenalty(costPt, edge_beta);
      currentPiece = currentPiece -> nxt;
    }
  }

  ///////////////////// CASE K != INF /////////////////////
  if(K != INFINITY)
  {
    ///Interval
    Interval new_interval = cost_intervalInterRoots(costPt, K);
    double AK = new_interval.geta();
    double BK = new_interval.getb();

    /// INTIALIZATION for Robust cost left and right
    coeff[0] = 0;
    coeff[1] = -a;
    coeff[2] = K + a * AK;
    Cost slopeLeftCost = Cost(coeff);  /// LEFT y = -ax + K + a * AK
    coeff[1] = a;
    coeff[2] = K  - a * BK;
    Cost slopeRightCost = Cost(coeff);  /// RIGHT y = ax + K - a * BK

    /// bounds of the tmp Piece
    double tmpA;
    double tmpB;
    int cas = 0;

    while(currentPiece != NULL)
    {
      tmpA = currentPiece -> m_interval.geta();
      tmpB = currentPiece -> m_interval.getb();

      if(tmpB <= AK){cas = 0;}
      if(BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;} // priority to AK over BK between tmpA and tmpB.

      switch(cas)
      {
        case 0 : currentPiece -> addCostAndPenalty(slopeLeftCost, edge_beta);
          break;
        case 1 : currentPiece -> addCostAndPenalty(slopeRightCost, edge_beta);
          break;
        case 2 : currentPiece -> addCostAndPenalty(costPt, edge_beta);
          break;
        case 3 :
        {
          // A) create nextPiece3 as a copy and update bound left
          Piece* nextPiece3 = new Piece(currentPiece);
          nextPiece3 -> m_interval.seta(BK); // changing interval bounds
          addCurrentPiecePlus1NotMove(nextPiece3);
          // B) update currentPiece cost and bound right
          currentPiece -> addCostAndPenalty(costPt, edge_beta); // adding costPt on the left
          currentPiece -> m_interval.setb(BK); // changing interval bounds
          break;
        }
        case 4 :
        {
          // A) create nextPiece4 as a copy and update bound left
          Piece* nextPiece4 = new Piece(currentPiece);
          nextPiece4 -> m_interval.seta(AK); // changing interval bounds
          addCurrentPiecePlus1NotMove(nextPiece4);
          // B) update currentPiece cost and bound right
          currentPiece -> addCostAndPenalty(slopeLeftCost, edge_beta); // adding slopeLeftCost on the left
          currentPiece -> m_interval.setb(AK); // changing interval bounds
          break;
        }
      }
      lastPiece = currentPiece; // DANGER: update lastPiece position
      move();
    }
  }
  delete[] coeff;
}



//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///
//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///

void ListPiece::LP_ts_Minimization(ListPiece& LP_edge)
{
  // Initialize LP_edge -> same range as this
  Interval newBounds = Interval(this -> head -> m_interval.geta(), this -> lastPiece -> m_interval.getb());
  LP_edge.setNewBounds(newBounds);

    //"MOST OF THE TIME" : Q2 > Q1
  Piece* Q1 = head;  /// first Piece to compare
  Piece* Q2 = LP_edge.head;  /// first Piece to compare

  Piece* Q12 = new Piece();
  Q12 -> m_interval = Interval(Q1 -> m_interval.geta(), Q1 -> m_interval.geta());
  int Bound_Q2_Minus_Q1 = 0;
  ///Q12 = Piece with an interval but no cost no label
  /// Bound_Q2_Minus_Q1
  /// = 0 if bound interval Q2 - bound interval Q1 == 0 : Q1 and Q2 stop
  /// = 1 if bound interval Q2 - bound interval Q1 > 0: Q2 stops
  /// = -1 if bound interval Q2 - bound interval Q1 < 0 : Q1 stops

  ///start info
  Piece* newHead = Q12;
  double M = lastPiece -> m_interval.getb(); //global right bound

  while(Q1 != NULL)
  {
    Bound_Q2_Minus_Q1 = -1;
    while(Bound_Q2_Minus_Q1 == -1)
    {
      /// right bound
      if(Q1 -> m_interval.getb() < Q2 -> m_interval.getb()){Bound_Q2_Minus_Q1 = 1;}
      if(Q1 -> m_interval.getb() == Q2 -> m_interval.getb()){Bound_Q2_Minus_Q1 = 0;}
      Q12 = Q12 -> pieceGenerator(Q1, Q2, Bound_Q2_Minus_Q1, M); ///add new Piece(s) to Q12
      if(Bound_Q2_Minus_Q1 < 1){Q2 = Q2 -> nxt;}
    }
    Q1 = Q1 -> nxt;
  }

  reset(); ///UPDATE ListPiece LP_ts[t + 1][i]
  head = newHead;
  currentPiece = newHead;
  lastPiece = Q12;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### operatorUp #####//////##### operatorUp #####//////##### operatorUp #####///
//##### operatorUp #####//////##### operatorUp #####//////##### operatorUp #####///

void ListPiece::operatorUp(ListPiece const& LP_state, unsigned int newLabel, unsigned int parentState)
{
  /// variable definition
  Piece* tmp; ///to follow LP_state List
  double currentValue; ///for the ListPiece to build, last current value
  double rightBound; ///value at the (right) bound of the last build interval. Up case
  bool constPiece; ///has the Piece to build constant cost?
  unsigned int counter = 1; ///number of the considered Piece in LP_edge
  Track trackUp = Track(newLabel, parentState, counter);
  Interval decreasingInterval; /// for interval building

  //////////////////
  ///First Piece head
  //////////////////
  head = new Piece();
  tmp = LP_state.head;

  /// INFO
  head -> m_info.setTrack(trackUp); ///set Track

  /// INTERVAL
  rightBound = tmp -> m_interval.geta(); ///tmp interval
  head -> m_interval.seta(rightBound);
  head -> m_interval.setb(rightBound);

  /// COST
  currentValue = cost_eval(tmp -> m_cost, rightBound); ///tmp cost value
  addConstant(head -> m_cost, currentValue);

  /// bool constPiece : is the first Piece constant? If cost increasing at bound, constPiece = true

  if(cost_argmin(tmp -> m_cost) <= rightBound && isConstant(tmp -> m_cost) == false){constPiece = true;}else{constPiece = false;}

  initializeCurrentPiece(); ///currentPiece = head
  ///////////////////////////

  while(tmp != NULL)
  {
    ///decreasingInterval for currentPiece to create based on current tmp
    decreasingInterval = tmp -> intervalMinLessUp(rightBound, currentValue, constPiece); ///"decreasing" interval
    decreasingInterval = decreasingInterval.intersection(tmp -> m_interval); ///decreasingInterval = intersection of decreasingInterval (=intervalMinLess) and interval of  tmp
    if(decreasingInterval.isEmpty() == false){trackUp.setPosition(counter);}

    /// paste new piece(s)
    currentPiece = currentPiece -> pastePieceUp(tmp, decreasingInterval, trackUp); ///add new Piece to BUILD
    ///

    ///UDPATES rightBound, currentValue, constPiece
    rightBound = currentPiece -> m_interval.getb(); ///new rightBound
    currentValue = cost_eval(currentPiece -> m_cost, rightBound); ///new currentValue (=the minimum)
    if(constPiece == true){if(decreasingInterval.isEmpty() == false){constPiece = false;}}
    if(constPiece == false){if(decreasingInterval.getb() < tmp -> m_interval.getb()){constPiece = true;}}

    tmp = tmp -> nxt;
    counter = counter + 1;
  }
  lastPiece = currentPiece;
}



//##### operatorDw #####//////##### operatorDw #####//////##### operatorDw #####///
//##### operatorDw #####//////##### operatorDw #####//////##### operatorDw #####///

void ListPiece::operatorDw(ListPiece const& LP_state, unsigned int newLabel, unsigned int parentState)
{
  /// variable definition
  Piece* tmp; ///to follow LP_state List
  double currentValue; ///for the ListPiece to build, last current value
  double leftBound; ///value at the (left) bound of the last build interval. Down case
  bool constPiece; ///has the Piece to build constant cost?
  unsigned int counter = 1; ///number of the considered Piece in LP_edge
  Track trackUp = Track(newLabel, parentState, counter);
  Interval decreasingInterval = Interval(); /// for interval building


  //////////////////
  ///First Piece head
  //////////////////
  head = new Piece();
  tmp = LP_state.head;

  /// INFO
  head -> m_info.setTrack(trackUp); ///set Track

  /// INTERVAL
  leftBound = tmp -> m_interval.getb();
  head -> m_interval.seta(leftBound);
  head -> m_interval.setb(leftBound);

  /// COST
  currentValue = cost_eval(tmp -> m_cost, leftBound);
  addConstant(head -> m_cost, currentValue);

  /// bool constPiece : is the first Piece constant? If cost increasing at bound, constPiece = true
  if(cost_argmin(tmp -> m_cost) >= leftBound && isConstant(tmp -> m_cost) == false){constPiece = true;}else{constPiece = false;}

  initializeCurrentPiece(); ///currentPiece = head

  ///////////////////////////

  while(tmp != NULL)
  {
    ///decreasingInterval for currentPiece to create based on current tmp
    decreasingInterval = tmp -> intervalMinLessDw(leftBound, currentValue, constPiece); ///"decreasing" interval
    //std::cout << leftBound << " " << currentValue << " " << constPiece << std::endl;
    //decreasingInterval.show();
    decreasingInterval = decreasingInterval.intersection(tmp -> m_interval); ///decreasingInterval = intersection of decreasingInterval (=intervalMinLess) and interval of  tmp
    if(decreasingInterval.isEmpty() == false){trackUp.setPosition(counter);}

    /// paste new piece(s)
    currentPiece = currentPiece ->pastePieceDw(tmp, decreasingInterval, trackUp); ///add new Piece to BUILD
    ///

    ///UDPATES rightBound, currentValue, constPiece
    leftBound = currentPiece -> m_interval.geta(); ///new rightBound
    currentValue = cost_eval(currentPiece -> m_cost, leftBound); ///new currentValue (=the minimum)
    if(constPiece == true){if(decreasingInterval.isEmpty() == false){constPiece = false;}}
    if(constPiece == false){if(decreasingInterval.geta() > tmp -> m_interval.geta()){constPiece = true;}}

    tmp = tmp -> nxt;
    counter = counter + 1;
  }
  lastPiece = currentPiece;
}



//####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######//
//####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######// //####### min_argmin_label_state_position_final #######//
///We test all the Piece

void ListPiece::get_min_argmin_label_state_position_ListPiece(double* response)
{
  Piece* tmp = head;

  ///INITIALIZATION
  tmp -> get_min_argmin_label_state_position(response);
  double current_min;

  tmp = tmp -> nxt;

  ///LOOP TESTS
  while(tmp != NULL)
  {
    current_min = cost_minInterval(tmp -> m_cost, tmp -> m_interval);
    if(current_min < response[0])
    {
      tmp -> get_min_argmin_label_state_position(response);
    }
    tmp = tmp -> nxt;
  }
}

//####### get_min_argmin_label_state_position_onePiece #######// //####### get_min_argmin_label_state_position_onePiece #######// //####### get_min_argmin_label_state_position_onePiece #######//
//####### get_min_argmin_label_state_position_onePiece #######// //####### get_min_argmin_label_state_position_onePiece #######// //####### get_min_argmin_label_state_position_onePiece #######//
///We test all the Piece

void ListPiece::get_min_argmin_label_state_position_onePiece(double* response, unsigned int position, Interval constrainedInterval, bool out, bool& forced)
{
  Piece* tmp = head;
  unsigned int nb = 1;
  while(nb < position){tmp = tmp -> nxt; nb = nb + 1;}
  tmp -> get_min_argmin_label_state_position(response);
  forced = false;

  /// argmin correction
  /// argmin correction
  /// argmin correction
  if(out == false)
  {
    if(constrainedInterval.isInside(response[1]) == false)
    {
      if(response[1] > constrainedInterval.getb()){response[1] = constrainedInterval.getb(); forced = true;}
      if(response[1] < constrainedInterval.geta()){response[1] = constrainedInterval.geta(); forced = true;}
    }
  }

  if(out == true)
  {
    if((constrainedInterval.geta() <= response[1]) && (response[1] <= constrainedInterval.getb()))
    {
      forced = true;
      if(response[1] - constrainedInterval.geta() < constrainedInterval.getb() - response[1]){response[1] = constrainedInterval.geta();}
      else{response[1] = constrainedInterval.getb();}
    }
  }
}

/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::operatorSum(ListPiece& LP1, ListPiece& LP2)
{
  //LP1 and LP2 must have the same bounds (min and max for parameter range)
  reset(); /// build a new ListPiece sum of LP1 and LP2
  LP1.initializeCurrentPiece();
  LP2.initializeCurrentPiece();

  ///BUILDING head pointer
  head = new Piece();
  currentPiece = head;
  currentPiece -> m_cost = addCost(LP1.currentPiece -> m_cost, LP2.currentPiece -> m_cost);
  currentPiece -> m_interval = LP1.currentPiece -> m_interval.intersection(LP2.currentPiece -> m_interval);
  //currentPiece -> m_info =

  if(LP1.currentPiece -> m_interval.getb() == LP2.currentPiece -> m_interval.getb()){LP1.move(); LP2.move();}
  if(LP1.currentPiece -> m_interval.getb() < LP2.currentPiece -> m_interval.getb()){LP1.move();}
  if(LP1.currentPiece -> m_interval.getb() > LP2.currentPiece -> m_interval.getb()){LP2.move();}

  while(LP1.currentPiece != NULL) // the  LP1 and LP2 currentPiece will be NULL at the same loop step
  {
    currentPiece -> nxt = new Piece();
    currentPiece = currentPiece -> nxt;

    currentPiece -> m_cost = addCost(LP1.currentPiece -> m_cost, LP2.currentPiece -> m_cost);
    currentPiece -> m_interval = LP1.currentPiece -> m_interval.intersection(LP2.currentPiece -> m_interval);
    //currentPiece -> m_info =

    if(LP1.currentPiece -> m_interval.getb() == LP2.currentPiece -> m_interval.getb()){LP1.move(); LP2.move();}
    if(LP1.currentPiece -> m_interval.getb() < LP2.currentPiece -> m_interval.getb()){LP1.move();}
    if(LP1.currentPiece -> m_interval.getb() > LP2.currentPiece -> m_interval.getb()){LP2.move();}

  }
  lastPiece = currentPiece;
}



//####### show #######// //####### show #######// //####### show #######//
//####### show #######// //####### show #######// //####### show #######//


void ListPiece::show() const
{
  //std::cout << "    HEAD      " << head << std::endl;
  //std::cout << "    CURRENTPI " << currentPiece << std::endl;
  //std::cout << "    LASTPIECE " << lastPiece << std::endl;
  head -> show();
}


//####### test #######// //####### test #######// //####### test #######//
//####### test #######// //####### test #######// //####### test #######//

void ListPiece::test()
{
  unsigned int t = 0;
  currentPiece = head;
  double leftEval;
  double rightEval;

  while(currentPiece != NULL)
  {
    t = t + 1;
    /// TEST INTERVAL
    if(currentPiece -> m_interval.getb() <= currentPiece -> m_interval.geta())
    {
      //std::cout << "DANGER DANGER DANGER DANGER DANGER DANGER" << std::endl;
      //std::cout << currentPiece -> m_interval.getb() << "  " << currentPiece -> m_interval.geta() << std::endl;
      //Rcpp::stop("Unexpected condition occurred: interval bounds");
    }

    /// TEST CONTINUITY
    if(currentPiece -> nxt != NULL)
    {
      rightEval = cost_eval(currentPiece -> m_cost, currentPiece -> m_interval.getb());
      leftEval = cost_eval(currentPiece -> nxt -> m_cost, currentPiece -> nxt -> m_interval.geta());
      if(fabs(rightEval - leftEval) > std::pow(10.0,-6.0))
      {
        //std::cout << std::endl;
        //std::cout << "          " << currentPiece;
        //std::cout << " #LABEL# "<< currentPiece -> m_info.getLabel() << " #STATE# " <<  currentPiece -> m_info.getState() << " POSITION " << currentPiece -> m_info.getPosition() << " ";
        //std::cout << " #INTERVAL# "<< currentPiece -> m_interval.geta() << " to " << currentPiece -> m_interval.getb() << " ";
        showCost(currentPiece -> m_cost);
        //std::cout << "          " << currentPiece -> nxt;
        //std::cout << " #LABEL# "<< currentPiece -> nxt -> m_info.getLabel() << " #STATE# " <<  currentPiece -> nxt -> m_info.getState() << " POSITION " << currentPiece -> nxt -> m_info.getPosition() << " ";
        //std::cout << " #INTERVAL# "<< currentPiece -> nxt -> m_interval.geta() << " to " << currentPiece -> nxt -> m_interval.getb() << " ";
        showCost(currentPiece -> nxt -> m_cost);
        //std::cout << " #rightEval# "<< rightEval << " #leftEval# " << leftEval << " #DIFF# " << leftEval - rightEval << std::endl;
        Rcpp::stop("Unexpected condition occurred: continuity error");
      }
    }

    currentPiece = currentPiece -> nxt;
  }

  //std::cout << "nb: " << t << std::endl;
}

