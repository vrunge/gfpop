#ifndef INTERVAL_H
#define INTERVAL_H

class Interval
{
  public:
    Interval();
    Interval(double a, double b);

    bool isEmpty() const;
    Interval intersection(Interval const& inter) const;

    void seta(double a);
    void setb(double b);
    double geta() const;
    double getb() const;

    bool isInside(double x) const;
    void axisSymmetry();

    double internPoint() const;

    void show() const;

  private:
    double m_a;
    double m_b;
};

#endif // INTERVAL_H
