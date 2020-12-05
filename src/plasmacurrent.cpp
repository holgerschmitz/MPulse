#include "plasmacurrent.hpp"

#include "../huerto/electromagnetics/fieldsolver.hpp"

#include <memory>

void PlasmaCurrentBlock::initParameters(schnek::BlockParameters &blockPars)
{
  CurrentBlock::initParameters(blockPars);

  blockPars.addParameter("em",&em,1.0);
  blockPars.addParameter("gamma",&gamma,0.01);
}

void PlasmaCurrentBlock::initCurrents(CurrentContainer &container)
{
  container.addCurrent(std::make_shared<PlasmaCurrent>(em, gamma, Z, gamma, boost::ref(*this)));
}

PlasmaCurrent::PlasmaCurrent(double em_, double mi_, double Z_, double gamma_, CurrentBlock &plasmaBlock_)
  : plasmaBlock(plasmaBlock_), em(em_), mi(mi_), Z(Z_), gamma(gamma_)
{}

void PlasmaCurrent::init()
{
  schnek::DomainSubdivision<Field> &subdivision = plasmaBlock.getContext().getSubdivision();

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  plasmaBlock.retrieveData("Ex", pEx);
  plasmaBlock.retrieveData("Ey", pEy);
  plasmaBlock.retrieveData("Ez", pEz);

  plasmaBlock.retrieveData("Rho", pRho);


  pJx = std::make_shared<Grid>(lowIn, highIn);
  pJy = std::make_shared<Grid>(lowIn, highIn);
  pJz = std::make_shared<Grid>(lowIn, highIn);

  plasmaBlock.addData("PlasmaJx", pJx);
  plasmaBlock.addData("PlasmaJy", pJy);
  plasmaBlock.addData("PlasmaJz", pJz);

}

void PlasmaCurrent::stepScheme(double dt)
{
  Field &Ex = *pEx;
  Field &Ey = *pEy;
  Field &Ez = *pEz;

  Field &Rho = *pRho;

  Grid &Jx = *pJx;
  Grid &Jy = *pJy;
  Grid &Jz = *pJz;

  Index low = Jx.getLo();
  Index high = Jx.getHi();

  const double gdtn = 1-0.5*gamma*dt;
  const double gdtd = 1+0.5*gamma*dt;
  const double emdt = dt*Z*em/mi;

#ifdef HUERTO_ONE_DIM
  for (int i=low[0]; i<high[0]; ++i) {
    double &jx = Jx(i);
    double &jy = Jy(i);
    double &jz = Jz(i);
    double rho = Rho(i);

    jx = (jx*gdtn - emdt*Ex(i)*rho)/gdtd;
    jy = (jy*gdtn - emdt*Ey(i)*rho)/gdtd;
    jz = (jz*gdtn - emdt*Ez(i)*rho)/gdtd;
  }
#endif

#ifdef HUERTO_TWO_DIM
  for (int i=low[0]; i<high[0]; ++i) {
    for (int j=low[1]; j<high[1]; ++j) {
      double &jx = Jx(i,j);
      double &jy = Jy(i,j);
      double &jz = Jz(i,j);
      double rho = Rho(i,j);

      jx = (jx*gdtn - emdt*Ex(i,j)*rho)/gdtd;
      jy = (jy*gdtn - emdt*Ey(i,j)*rho)/gdtd;
      jz = (jz*gdtn - emdt*Ez(i,j)*rho)/gdtd;
    }
  }
#endif

#ifdef HUERTO_THREE_DIM
  for (int i=low[0]; i<high[0]; ++i) {
    for (int j=low[1]; j<high[1]; ++j) {
      for (int k=low[2]; k<high[2]; ++k) {
        double &jx = Jx(i,j,k);
        double &jy = Jy(i,j,k);
        double &jz = Jz(i,j,k);
        double rho = Rho(i,j,k);

        jx = (jx*gdtn - emdt*Ex(i,j,k)*rho)/gdtd;
        jy = (jy*gdtn - emdt*Ey(i,j,k)*rho)/gdtd;
        jz = (jz*gdtn - emdt*Ez(i,j,k)*rho)/gdtd;
      }
    }
  }
#endif

}

