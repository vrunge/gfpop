#include "Track.h"

Track::Track()
{
  myLabel = 0;
  myParentState = 0;
  myParentPosition = 0;
}

Track::Track(unsigned int label, unsigned int state, unsigned int position)
{
  myLabel = label;
  myParentState = state;
  myParentPosition = position;
}

unsigned int Track::getLabel() const {return(myLabel);}
unsigned int Track::getState() const {return(myParentState);}
unsigned int Track::getPosition() const {return(myParentPosition);}

void Track::setPosition(unsigned int position){myParentPosition = position;}
void Track::setTrack(Track const& newTrack)
{
  myLabel = newTrack.getLabel();
  myParentState = newTrack.getState();
  myParentPosition = newTrack.getPosition();
}

void Track::reversePosition(unsigned int length){myParentPosition = length - myParentPosition + 1;}
