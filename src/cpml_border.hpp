#ifndef MPULSE_CPML_BORDER_H
#define MPULSE_CPML_BORDER_H

#include "mpulse.hpp"
#include "../huerto/electromagnetics/current.hpp"

class CPMLBorder : public CurrentBlock
{
  public:
    void initCurrents(CurrentContainer &container);
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void init();
  private:

    void initCoefficients();

    int thickness;
    double kappaMax;
    double aMax;
    double sigmaMax;
    double eps;
};

class CPMLBorderCurrent : public Current
{
  public:
    CPMLBorderCurrent(int thickness, Direction dir, bool isH,
                      double kappaMax, double aMax, double sigmaMax, double eps,
                      CurrentBlock &borderBlock);
  protected:
    bool reverse;
    int thickness;

    int dim;
    int transverse1, transverse2;

    Direction dir;
    bool isH;
    int lowOffset;
    int highOffset;

    int zerolayer;

    double kappaMax;
    double aMax;
    double sigmaMax;
    double eps;

    Grid1d bCoeff;
    Grid1d cCoeff;

    CurrentBlock &borderBlock;

    void makeCoeff();
};

class CPMLBorderECurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderECurrent( int thickness_, Direction dir_,
                        double kappaMax_, double aMax_, double sigmaMax_, double eps_,
                        CurrentBlock &borderBlock_);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:

    pField pB[3];
    pGrid pPsi[2];
    double dx;
};

class CPMLBorderHCurrent : public CPMLBorderCurrent
{
  public:
    CPMLBorderHCurrent( int thickness_, Direction dir_,
                        double kappaMax_, double aMax_, double sigmaMax_, double eps_,
                        CurrentBlock &borderBlock_);

    void init();

    void stepSchemeInit(double dt);
    void stepScheme(double dt);
  protected:
    pField pE[3];
    pGrid pPsi[2];
    double dx;
};


#endif
