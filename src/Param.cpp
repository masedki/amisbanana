#include "Param.h"

Param::Param(){
  m_pi = ones<vec>(0);
  m_mu = ones<mat>(0,0);
  m_s = ones<cube>(0,0,0);
}

Param::Param(const Param & param){
  m_pi = param.m_pi;  
  m_mu = param.m_mu;
  m_s = param.m_s;
}

Param::Param(const Data * x_p, const int & g){
  m_mu = ones<mat>(x_p->m_x.n_cols, g);
  m_s = zeros<cube>(x_p->m_x.n_cols, x_p->m_x.n_cols, g);
  m_s.each_slice() += eye<mat>(x_p->m_x.n_cols, x_p->m_x.n_cols);
  m_pi = ones<vec>(g)/g;  
  for(int k = 0; k < g; ++k){
      ivec who = randi<ivec>(1, distr_param(0, x_p->m_x.n_rows -1));
      m_mu.col(k) = trans(x_p->m_x.row(who(0)));
     }
}

