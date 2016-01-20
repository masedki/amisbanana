#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef Data_H
#define Data_H

class Data{
  public:
  Mat<double> m_x;  
  
  Data(){};
  Data(const Mat<double> &);
  ~Data(){};
  
};
#endif