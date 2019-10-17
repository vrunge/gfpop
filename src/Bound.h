#ifndef BOUND_H
#define BOUND_H


class Bound
{
  public:
    Bound();
    Bound(double mini, double maxi, bool isConstr = false);

    double getm() const;
    double getM() const;
    bool getIsConstrained() const;

    void setm(double mini);
    void setM(double maxi);
    void setIsConstrained(bool isConstr);

  private:
    double m; ///minimum value of the theta interval
    double M; ///maximum value of the theta interval
    bool isConstrained; ///if false, all the data are inside the interval [m,M]


};

#endif // BOUND_H
