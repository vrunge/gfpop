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
  while(head != NULL)
  {
    Piece* pieceToDelete = head;
    head = head -> nxt;
    delete(pieceToDelete);
  }
}

//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///
//##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####//////##### setUniquePieceCostToInfinity #####///

void ListPiece::setUniquePieceCostToInfinity()
{
  head -> m_cost.constant = INFINITY;
}

//##### reset #####//////##### reset #####//////##### reset #####///
//##### reset #####//////##### reset #####//////##### reset #####///

void ListPiece::reset()
{
  while(head != NULL)
  {
    Piece* pieceToDelete = head;
    head = head -> nxt;
    delete(pieceToDelete);
  }
  currentPiece = NULL;
  lastPiece = NULL;
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
    length = length + 1;
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


//##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####///
//##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####//////##### addCurrentPiecePlus1 #####///

void ListPiece::addCurrentPiecePlus1(Piece* newPiece)
{
  newPiece -> nxt = currentPiece -> nxt;
  currentPiece -> nxt = newPiece;
}

//##### addNewLastPiece #####//////##### addNewLastPiece #####//////##### addNewLastPiece #####///
//##### addNewLastPiece #####//////##### addNewLastPiece #####//////##### addNewLastPiece #####///

void ListPiece::addNewLastPiece(Piece* newPiece)
{
  lastPiece -> nxt = newPiece;
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

//##### copy #####//////##### copy #####//////##### copy #####///
//##### copy #####//////##### copy #####//////##### copy #####///

void ListPiece::copy(ListPiece const& LP_edge)
{
  Piece* tmp = LP_edge.head;
  head = tmp;

  while(tmp != NULL)
  {
    currentPiece = tmp -> copy(); //copy content in Piece
    currentPiece -> nxt = tmp -> nxt; //copy nxt pointer
    tmp = tmp -> nxt;
  }
  lastPiece = currentPiece;
}


//##### shift #####//////##### shift #####//////##### shift #####///
//##### shift #####//////##### shift #####//////##### shift #####///

void ListPiece::shift(double parameter)
{
  currentPiece = head;
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    Interval inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interShift(inter.geta(), parameter));
    currentPiece -> m_interval.setb(cost_interShift(inter.getb(), parameter));
    ///MOVE Cost
    cost_shift(currentPiece -> m_cost, parameter);
    currentPiece = currentPiece -> nxt;
  }
}


//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///
//##### expDecay #####//////##### expDecay #####//////##### expDecay #####///

void ListPiece::expDecay(double gamma)
{
  currentPiece = head;
  while(currentPiece != NULL)
  {
    ///MOVE bounds
    Interval inter = currentPiece -> m_interval;
    currentPiece -> m_interval.seta(cost_interExpDecay(inter.geta(), gamma));
    currentPiece -> m_interval.setb(cost_interExpDecay(inter.getb(), gamma));
    ///MOVE Cost
    cost_expDecay(currentPiece -> m_cost, gamma);
    currentPiece = currentPiece -> nxt;
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///
//##### LP_edges_constraint #####//////##### LP_edges_constraint #####//////##### LP_edges_constraint #####///

void ListPiece::LP_edges_constraint(ListPiece const& LP_state, Edge const& edge, unsigned int newLabel)
{
  /// build a new LP_edges from scratch
  reset();

  /// 4 types of edges : null, std, up, down
  ///
  /// EDGE PARAMETERS
  ///
  std::string edge_ctt = edge.getConstraint();
  double edge_parameter = edge.getParameter(); /// always positive
  double edge_beta = edge.getBeta();
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
    double mini = INFINITY;
    double getmin;
    unsigned int counter = 0;
    unsigned int counterMini;

    ///find the minimum
    while(tmp != NULL)
    {
      counter = counter + 1;
      getmin = cost_minInterval(tmp -> m_cost, tmp -> m_interval);
      if(getmin < mini){mini = getmin; counterMini = counter;}
      tmp = tmp -> nxt;
    }

    ///add onePiece to LP_edges
    Piece* onePiece = new Piece();
    onePiece ->addCostAndPenalty(Cost(), mini + edge_beta); /// Cost() = 0
    onePiece -> m_info = Track(newLabel, parentState, counterMini);
    addNewLastPiece(onePiece);

  }

  //################
  if(edge_ctt == "up")
  {
    operatorUpDown(LP_state, newLabel, parentState, true);
    if(edge_parameter > 0){shift(edge_parameter);} ///edge_parameter = left/right decay
  }

  //################
  if(edge_ctt == "down")
  {
    ///LP_state working copy
    ListPiece LP_stateCopy = ListPiece();
    LP_stateCopy.copy(LP_state);
    ///reverse LP_stateCopy
    unsigned int length = 0;
    LP_stateCopy.reverseAndCount(length);
    ///down operations and reverse result
    operatorUpDown(LP_stateCopy, newLabel, parentState, false);
    reverseAndSetTrackPosition(length);

    if(edge_parameter > 0){shift(-edge_parameter);} ///edge_parameter = left/right decay
  }

}

//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///
//##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####//////##### LP_edges_addPointAndPenalty #####///

void ListPiece::LP_edges_addPointAndPenalty(Edge const& edge, Point const& pt)
{
  /// get edge data ///
  double K = edge.getKK();
  double a = edge.getAA();
  double penalty = edge.getBeta();

  initializeCurrentPiece();

  ///////////////////// CASE K == INF /////////////////////
  ///////////////////// CASE K == INF /////////////////////
  if(K == INFINITY)
  {
    Cost costPt = Cost(cost_coeff(pt));
    while(currentPiece != NULL)
    {
      currentPiece -> addCostAndPenalty(costPt, penalty);
      move();
    }
  }

  ///////////////////// CASE K != INF /////////////////////
  ///////////////////// CASE K != INF /////////////////////
  if(K != INFINITY)
  {
    /// INTIALIZATION
    ///Interval
    Cost costInter = Cost(cost_coeff(pt));
    Interval new_interval = cost_intervalInterRoots(costInter, K);

    /// Putting the bounds in variables AK and BK
    double AK = new_interval.geta();
    double BK = new_interval.getb();

    double* coeff = new double[3];
    coeff[0] = 0;
    coeff[1] = -a;
    coeff[2] = K + a * AK;
    Cost slopeLeftCost = Cost(coeff);  /// LEFT y = -ax + K + a * AK
    coeff[1] = a;
    coeff[2] = K  - a * BK;
    Cost slopeRightCost = Cost(coeff);  /// RIGHT y = ax + K - a * BK
    Cost costPt = Cost(cost_coeff(pt));  /// CENTER pt

    /// bounds of the tmp Piece
    double tmpA = 0;
    double tmpB = 0;

    int cas = 0;

    while(currentPiece != NULL)
    {
      tmpA = currentPiece -> m_interval.geta();
      tmpB = currentPiece -> m_interval.getb();

      if(tmpB <= AK){cas = 0;}
      if(BK <= tmpA){cas = 1;}
      if(AK <= tmpA && tmpB <= BK){cas = 2;}
      if(tmpA < BK && BK < tmpB){cas = 3;}
      if(tmpA < AK && AK < tmpB){cas = 4;}///priority to AK over BK between tempA and tempB.

      switch(cas)
      {
      case 0 : currentPiece -> addCostAndPenalty(slopeLeftCost, penalty);
        break;
      case 1 : currentPiece -> addCostAndPenalty(slopeRightCost, penalty);
        break;
      case 2 : currentPiece -> addCostAndPenalty(costPt, penalty);
        break;
      case 3 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = currentPiece -> copy();
        addCurrentPiecePlus1(new_piece);
        ///adding costPt on the left
        currentPiece -> addCostAndPenalty(costPt, penalty);
        ///changing interval bounds
        currentPiece -> m_interval.setb(BK);
        new_piece -> m_interval.seta(BK);
        break;
      }

      case 4 :
      {
        ///copying currentPiece in new_piece add adding new_piece after currentPiece
        Piece* new_piece = new Piece(currentPiece);
        addCurrentPiecePlus1(new_piece);
        ///adding slopeLeftCost on the left
        currentPiece -> addCostAndPenalty(slopeLeftCost, penalty);
        ///changing interval bounds
        currentPiece -> m_interval.setb(AK);
        new_piece -> m_interval.seta(AK);
        break;
      }
      }
    move();
    }
  }


}



//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///
//##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####//////##### LP_ts_Minimization #####///

void ListPiece::LP_ts_Minimization(ListPiece const& LP_edge)
{


}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//##### operatorUpDown #####//////##### operatorUpDown #####//////##### operatorUpDown #####///
//##### operatorUpDown #####//////##### operatorUpDown #####//////##### operatorUpDown #####///

void ListPiece::operatorUpDown(ListPiece const& LP_state, unsigned int newLabel, unsigned int parentState, bool upDirection)
{
  /// variable definition
  Piece* tmp; ///to follow LP_edge List
  double currentValue; ///for the ListPiece to build, last current value
  double rightBound; ///value at the (right) bound of the last build interval. Up case
  double leftBound; ///value at the (left) bound of the last build interval. Down case
  bool constPiece; ///has the Piece to build constant cost?
  unsigned int counter = 1; ///number of the considered Piece in LP_edge
  Track trackUp = Track(newLabel, parentState, counter);
  Interval decreasingInterval = Interval(); /// for interval building

  ///upDirection///upDirection///upDirection///upDirection///upDirection///upDirection///upDirection
  if(upDirection == true)
  {
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
    addmyConstant(head -> m_cost, currentValue);

    /// bool constPiece : is the first Piece constant? If cost increasing at bound, constPiece = true
    if(cost_argmin(tmp -> m_cost) <= rightBound){constPiece = true;}else{constPiece = false;}

    initializeCurrentPiece(); ///currentPiece = head

    ///////////////////////////

    while(tmp != NULL)
    {
      ///decreasingInterval for currentPiece to create based on current tmp
      decreasingInterval = tmp -> intervalMinLess(rightBound, currentValue, constPiece, true); ///"decreasing" interval
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
  }

  ///downDirection///downDirection///downDirection///downDirection///downDirection///downDirection///downDirection
  if(upDirection == false)
  {
    //////////////////
    ///First Piece head
    //////////////////
    head = new Piece();
    tmp = LP_state.head;

    /// INFO
    head -> m_info.setTrack(trackUp); ///set Track

    /// INTERVAL
    leftBound = tmp -> m_interval.getb();
    head -> m_interval.seta(rightBound);
    head -> m_interval.setb(rightBound);

    /// COST
    currentValue = cost_eval(tmp -> m_cost, leftBound);
    addmyConstant(head -> m_cost, currentValue);

    /// bool constPiece : is the first Piece constant? If cost increasing at bound, constPiece = true
    if(cost_argmin(tmp -> m_cost) <= leftBound){constPiece = true;}else{constPiece = false;}

    initializeCurrentPiece(); ///currentPiece = head

    ///////////////////////////

    while(tmp != NULL)
    {
      ///decreasingInterval for currentPiece to create based on current tmp
      decreasingInterval = tmp -> intervalMinLess(leftBound, currentValue, constPiece, false); ///"decreasing" interval
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
  }

}


/////////////////////////////////////////
/////////////////////////////////////////


void ListPiece::show()
{
  initializeCurrentPiece();

  while(currentPiece != NULL)
  {
    currentPiece -> show();
    move();
  }
}

