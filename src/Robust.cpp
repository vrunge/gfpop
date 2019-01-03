#include "Robust.h"

#include<iostream>

Robust::Robust(double K, double a) : threshold(K), slope(a){}

double Robust::getThreshold() const {return(threshold);}
double Robust::getSlope() const {return(slope);}
std::string Robust::getRobustType() const {return(robustType);}


void Robust::setThreshold(double K){threshold = K;}
void Robust::setSlope(double a){slope = a;}


// ### findRobustType ### ////// ### findRobustType ### ////// ### findRobustType ### ///
// ### findRobustType ### ////// ### findRobustType ### ////// ### findRobustType ### ///

void Robust::findRobustType()
{
  if(threshold == INFINITY){robustType = "notRobust";}

  if(threshold != INFINITY)
  {
    if(slope == 0){robustType = "biweight";}
    if(slope > 0){robustType = "Huber";}
  }
}

// ### show ### ////// ### show ### ////// ### show ### ////// ### show ### ///
// ### show ### ////// ### show ### ////// ### show ### ////// ### show ### ///

void Robust::show() const
{
  //std::cout<< "ROBUST" << " -- robustType : " << robustType << std::endl;
  //std::cout<< "K : " << threshold << std::endl;
  //std::cout<< "a : " << slope << std::endl;
}

