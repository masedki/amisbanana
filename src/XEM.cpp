#include "XEM.h"
const double log2pi = std::log(2.0 * M_PI);


Col<double> ldtmvt(const mat & x, const Col<double> & mu, const mat & s_sig, const int & nu){
  mat s_chol  = zeros<mat>(x.n_cols, x.n_cols);
  Col<double> ldensity(x.n_rows);
  const bool success =  chol(s_chol, s_sig);
  if(success==false)
  {
    Col<double> ldensity=log(0.)*ones<vec>(x.n_rows);
  }
  else{
    const double slogdet = 2 * sum(log(s_chol.diag())), d = s_sig.n_cols,  myconst = Rf_lgammafn((nu+d)/2.)  - Rf_lgammafn(nu/2.) - d*log(M_PI * nu)/2. - slogdet/2. ;
    mat s_inv = inv_sympd(s_sig);

    for(int i = 0; i < x.n_rows; ++i)
      ldensity(i) = myconst - 0.5*(nu+d)*log(1+ as_scalar((x.row(i)-trans(mu))*(s_inv*(trans(x.row(i))-mu)))/nu);
  }
  return(ldensity);
};



XEM::XEM(const S4 * popEM, const S4 * amisstrategy)
{

  x_p = new Data(as<mat>(popEM->slot("x")));
  g = as<int>(amisstrategy->slot("g"));
  nu = as<int>(amisstrategy->slot("nu"));
  nbSmall= as<int>(amisstrategy->slot("nbSmall"));
  iterSmall = as<int>(amisstrategy->slot("iterSmall"));
  nbKeep = as<int>(amisstrategy->slot("nbKeep"));
  iterKeep = as<int>(amisstrategy->slot("iterKeep"));
  tolKeep = as<double>(amisstrategy->slot("tolKeep"));
  iterCurrent= iterSmall;
  m_nbdegenere=0;
  tmplogproba = zeros<mat>(x_p->m_x.n_rows, g);
  maxtmplogproba = zeros<vec>(x_p->m_x.n_rows);
  loglikSmall = ones<vec>(nbSmall)*log(0.);
  loglikoutput = log(0.);
  for(int i=0;i<nbSmall;i++)
    paramCand.push_back(Param(x_p, g));
  rowsums = ones<vec>(x_p->m_x.n_rows);
}

int XEM::FiltreDegenere(){
  int output = 0;
  bool success = false;
  Col<double> s_det(g);
  Mat<double> s_chol(x_p->m_x.n_cols, x_p->m_x.n_cols);
  for(int k = 0; k<g; ++k)
  {
    bool success =  chol(s_chol, paramCurrent_p->m_s.slice(k));
    if(success==false)
      s_det(k) = 0.;
    else
      s_det(k) =  exp(2 * sum(log(s_chol.diag())));

  }
  if (min(s_det)<0.00001)
    output = 1;
  return output;
}


void XEM::ComputeTmpLogProba()
{
  tmplogproba = zeros<mat>(x_p->m_x.n_rows, g);
  for(int k=0; k<g; ++k)
    tmplogproba.col(k) += zeros<vec>(x_p->m_x.n_rows) + log(paramCurrent_p->m_pi(k)) + ldtmvt(x_p->m_x, paramCurrent_p->m_mu.col(k), paramCurrent_p->m_s.slice(k), nu);

}
void XEM::ComputeTmpU()
{
  mat s_inv =  zeros<mat>(x_p->m_x.n_cols, x_p->m_x.n_cols);
  //mat s_k = zeros<mat>(x_p->m_x.n_cols, x_p->m_x.n_cols) ;
  tmpu = zeros<mat>(x_p->m_x.n_rows, g);
  for(int k=0; k<g; ++k)
  {
    const bool success = inv_sympd(s_inv, paramCurrent_p->m_s.slice(k));
    if(success!=false)
//     while(success == false)
//     {
//       success = inv_sympd(s_inv, s_k);
//       if(success == false)
//       {
//         s_k += eye(s_k.n_rows,s_k.n_rows) * 1e-6;
//       }
//     }

    for(int i=0; i<x_p->m_x.n_rows; ++i)
      tmpu(i,k) = (nu+x_p->m_x.n_cols)/(nu+as_scalar((x_p->m_x.row(i)-trans(paramCurrent_p->m_mu.col(k)))*(s_inv*(trans(x_p->m_x.row(i))-paramCurrent_p->m_mu.col(k)))));
  }
}

double XEM::ComputeLogLik(){
  ComputeTmpLogProba();
  maxtmplogproba = max(tmplogproba, 1);
  double output=0;
  if (min(maxtmplogproba) == 0){
    output = log(0.);
  }else{
    for (int k=0; k<g; k++) tmplogproba.col(k)-=maxtmplogproba;
    tmplogproba = exp(tmplogproba);
    rowsums = sum(tmplogproba,1);
    output = sum(maxtmplogproba) + sum(log(rowsums));
  }
  //cout<<"loglik ---->"<< output <<endl;
  return output;
}


void XEM::Estep(){
  ComputeTmpU();
  for(int k=0; k<g; k++) tmplogproba.col(k) = tmplogproba.col(k)/rowsums;
}
void XEM::OneEM(){
  double loglike = ComputeLogLik(), prec = log(0.);
  int it=0;
  while ( (it<iterCurrent) && ((loglike-prec)>tolKeep) ){
    it ++;
    Estep();
    Mstep();
    prec = loglike;
    loglike = ComputeLogLik();
  }

}
void XEM::Mstep()
{
  mat tmpprobau = tmplogproba % tmpu;
  paramCurrent_p->m_pi = trans(mean(tmplogproba, 0));
  for(int k = 1; k < g; ++k)
    for(int j =1; j<x_p->m_x.n_cols; ++j)
      paramCurrent_p->m_mu(j,k) = sum(tmpprobau.col(k) % x_p->m_x.col(j))/sum(tmpprobau.col(k));
  paramCurrent_p->m_s = zeros<cube>(x_p->m_x.n_cols, x_p->m_x.n_cols, g);

  for(int k=0; k<g; ++k)
  {
    for(int i=0; i<x_p->m_x.n_rows; ++i)
      paramCurrent_p->m_s.slice(k)+= tmpprobau(i,k)*((trans(x_p->m_x.row(i)) - paramCurrent_p->m_mu.col(k))*(x_p->m_x.row(i)-trans(paramCurrent_p->m_mu.col(k))));
    paramCurrent_p->m_s.slice(k)*=(1/sum(tmplogproba.col(k)));
  }

}

void XEM::SwitchParamCurrent(int ini){paramCurrent_p = &paramCand[ini];}

void XEM::Run(){
  // Partie Small EM
  int degenere = 0;
  for (int ini=0; ini<nbSmall; ini++){
    //cout << "ini -->"<< ini <<endl;
    SwitchParamCurrent(ini);
    OneEM();
    loglikSmall(ini) = ComputeLogLik();
    if (loglikSmall(ini) != loglikSmall(ini))
      loglikSmall(ini) = -999999999999;
  }
  // On conserve les meilleurs initialisations
  uvec indices = sort_index(loglikSmall);
  iterCurrent = iterKeep;
  m_nbdegenere = 0;
  //if (nbSmall > nbKeep)    loglikeSmall( indices.head(nbSmall - nbKeep) ) = loglikeSmall( indices.head(nbSmall - nbKeep) ) + log(0);

  for (int tmp1=0; tmp1<nbKeep; tmp1++){
    //cout << "tmp1 -->"<< tmp1 <<endl;
    SwitchParamCurrent(indices(nbSmall - tmp1 - 1));
    OneEM();
    loglikSmall(indices(nbSmall - tmp1 - 1)) = ComputeLogLik();
    if (loglikSmall(indices(nbSmall - tmp1 - 1)) != loglikSmall(indices(nbSmall - tmp1 - 1))){
      m_nbdegenere ++;
      loglikSmall(indices(nbSmall - tmp1 - 1)) = -999999999999;
    }
  }
  uword  index;
  double indicebest = (loglikSmall).max(index);
  SwitchParamCurrent(index);
  loglikoutput = ComputeLogLik();
  degenere = FiltreDegenere();
  if (degenere==1){
    m_nbdegenere ++;
    loglikoutput = -999999999999;
  }
}

void XEM::Output(S4 * reference_p){
  if (m_nbdegenere<nbKeep){
    as<S4>(reference_p->slot("param")).slot("mu") = paramCurrent_p->m_mu;
    as<S4>(reference_p->slot("param")).slot("pi") = paramCurrent_p->m_pi;
    as<S4>(reference_p->slot("param")).slot("nu") = nu;
    as<S4>(reference_p->slot("param")).slot("sigma") = wrap(paramCurrent_p->m_s);
    as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") =  double(m_nbdegenere)/double(nbKeep);
    as<S4>(reference_p->slot("criteria")).slot("loglik") =  loglikoutput;
  }
  else{
    as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") = 1;
  }


}


// [[Rcpp::export]]
void runXEM(S4  reference) {
  S4 * reference_p = &reference;
  S4  popEM =  as<S4>(reference.slot("pop")), amisstrategy = as<S4>(reference.slot("strategy"));
  XEM *xem_p  = new XEM(&popEM, &amisstrategy);
  xem_p->Run();
  xem_p->Output(reference_p);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::Mat<double> rbanana(const int n, const int d, const double b)
{

  mat Y = arma::randn(n, d);
  Y.col(0) = 10.0*Y.col(0);
  Y.col(1)  = Y.col(1) -  b*(pow(Y.col(0),2) -  100);
  return Y;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::Col<double> dbanana(arma::Mat<double> & x, const double & b){
  mat clonex (x);
  clonex.col(1) += b*pow(clonex.col(0),2)- b*100;
  clonex.col(0)  = clonex.col(0)/10.;
  clonex = -0.5*(log2pi + pow(clonex, 2));
  return(0.1*exp(sum(clonex,1)));

}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void computedensities(S4 reference, const double b)
{
  S4 * reference_p = &reference;
  Mat<double> tmp_mu  = as<mat>(as<S4>(reference_p->slot("param")).slot("mu"));
  Mat<double> tmp_pi  = as<mat>(as<S4>(reference_p->slot("param")).slot("pi"));
  Cube<double> tmp_s  = as<cube>(as<S4>(reference_p->slot("param")).slot("sigma"));
  mat x = as<mat>(as<S4>(reference_p->slot("pop")).slot("x"));
  int nu = as<int>(as<S4>(reference_p->slot("param")).slot("nu"));
  Col<double> tmpd = zeros<vec>(x.n_rows);
  for(int k=0; k<tmp_mu.n_cols; ++k)
    tmpd += exp(zeros<vec>(x.n_rows) + log(tmp_pi(k)) + ldtmvt(x, tmp_mu.col(k), tmp_s.slice(k),  nu));

  Col<double> tmpdb = dbanana(x, b);
  as<S4>(reference_p->slot("pop")).slot("denom") = wrap(tmpd);
  as<S4>(reference_p->slot("pop")).slot("numera") = wrap(tmpdb);
  as<S4>(reference_p->slot("pop")).slot("w") = wrap(tmpdb/tmpd);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void computedmixture(S4 pop, S4 param, int n)
{
  S4 * pop_p = &pop;
  S4 * param_p = &param;
  Mat<double> tmp_mu  = as<mat>(param_p->slot("mu"));
  Mat<double> tmp_pi  = as<mat>(param_p->slot("pi"));
  Cube<double> tmp_s  = as<cube>(param_p->slot("sigma"));
  mat x = as<mat>(pop_p->slot("x"));
  int nu = as<int>(param_p->slot("nu"));
  Col<double> tmpd(as<vec>(pop_p->slot("denom")));
  for(int k=0; k<tmp_mu.n_cols; ++k)
    tmpd += n*exp(zeros<vec>(x.n_rows) + log(tmp_pi(k)) + ldtmvt(x, tmp_mu.col(k), tmp_s.slice(k),  nu));
  pop_p->slot("denom") = wrap(tmpd);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::Col<double> dlogistic(arma::Mat<double> & x, arma::Col<double> & s)
{
  mat clonex(x);

  for(int j=0; j<clonex.n_cols; ++j)
  {
    clonex.col(j) = exp(-clonex.col(j)/s(j));
    clonex.col(j) = clonex.col(j)/(s(j)*pow(1+clonex.col(j), 2));
  }
  return prod(clonex,1);
}
