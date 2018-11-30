#include "Track.h"

Track::Track()
{
  myLabel = 0;
  myParentState = -1;
  myParentPosition = -1;
}

Track::Track(int label, int state, int position)
{
  myLabel = label;
  myParentState = state;
  myParentPosition = position;
}


int Track::getLabel() const {return(myLabel);}
int Track::getState() const {return(myParentState);}
int Track::getPosition() const {return(myParentPosition);}

void Track::setLabel(int label){myLabel = label;}
void Track::setState(int state){myParentState = state;}
void Track::setPosition(int position){myParentPosition = position;}

void Track::setTrack(int label, int state, int position)
{
  myLabel = label;
  myParentState = state;
  myParentPosition = position;
}


void Track::setTrack(Track const& newTrack)
{
  myLabel = newTrack.getLabel();
  myParentState = newTrack.getState();
  myParentPosition = newTrack.getPosition();
}



void Track::axisSymmetry(int length){myParentPosition = length - myParentPosition + 1;}
