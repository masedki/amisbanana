// #include <iostream>
// #include <iomanip>
// #include <RcppArmadillo.h>
// #include <Rcpp.h>
// 
// using namespace std;
// using namespace arma;
// using namespace Rcpp;




#ifndef Param_H
#define Param_H
#include "Data.h"

class Param{
  public:
  Col<double> m_pi;
  Mat<double> m_mu;
  Cube<double> m_s;
  
  
  Param();
  Param(const Param & param);
  Param(const Data *, const int &);
  ~Param(){};
};
#endif