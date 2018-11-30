#ifndef TRACK_H
#define TRACK_H

class Track
{
  public:
    Track();
    Track(int label, int state, int position);

    int getLabel() const;
    int getState() const;
    int getPosition() const;

    void setLabel(int label);
    void setState(int state);
    void setPosition(int position);
    void setTrack(int label, int state, int position);
    void setTrack(Track const& newTrack);

    void axisSymmetry(int length); ///start counting the label from the end of Piece

  private:
    int myLabel; ///label of the Piece
    int myParentState; ///parent state of the Piece
    int myParentPosition; ///position of the parent Piece in the list of Pieces in functional cost

};

#endif // TRACK_H
