//
//  XEM.h
//
//
//  Created by sedki on 15/01/2016.
//
//

#ifndef _XEM_h
#define _XEM_h
#include "Data.h"
#include "Param.h"
//Col<double> ldtmvt(mat &, vec &, mat &, const int &)

class XEM{
public:
    int g, nu, nbSmall, iterSmall, nbKeep, iterKeep, iterCurrent, m_nbdegenere;
    Col<double> maxtmplogproba, loglikSmall,  rowsums;
    Mat<double> tmplogproba, tmpu;
    vector<Param> paramCand;
    Param * paramCurrent_p;
    const Data  * x_p;
    double tolKeep, loglikoutput;
    XEM(const S4 * , const S4 *); // Constructor
    XEM(); // Constructor
    ~XEM(){} ; // Destructor
    void Output(S4 * reference_p); // gestion des sorties
    double loglik();
    void ComputeTmpLogProba();
    void ComputeTmpU();
    double ComputeLogLik();
    void Estep();
    void Mstep();
    int  FiltreDegenere();
    void SwitchParamCurrent(int);
    void OneEM();
    void Run();
    //void smallEM(void);
    //void gettik(void);
    //void Gettau(void);
    //void Updatemu(void);
    //void GetEmpiricalCovariance(void);
    //void Updatesigma(void);
    //int GetNbClust(void);
    
};
#endif