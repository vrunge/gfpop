#ifndef TRACK_H
#define TRACK_H

class Track
{
  public:
    Track();
    Track(unsigned int label, unsigned int state, unsigned int position);

    unsigned int getLabel() const;
    unsigned int getState() const;
    unsigned int getPosition() const;

    void setLabel(unsigned int label);
    void setState(unsigned int state);
    void setPosition(unsigned int position);
    void setTrack(unsigned int label, unsigned int state, unsigned int position);
    void setTrack(Track const& newTrack);

  private:
    unsigned int myLabel; ///label of the Piece
    unsigned int myParentState; ///parent state of the Piece
    unsigned int myParentPosition; ///position of the parent Piece in the list of Pieces in functional cost

};

#endif // TRACK_H
