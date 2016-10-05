/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/grid.hpp>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

static const size_t DIMENSION = 3;

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Field<double, DIMENSION> Grid;
typedef boost::shared_ptr<Grid> pGrid;
typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );

static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

class MPulse : public schnek::Block, schnek::BlockContainer<FieldSolver>
{
  private:
    static Index globalMax;
    Index gridSize;
    Vector size;
    Vector dx;
    double dt;

    int tMax;
    Range innerRange;
    pGrid Ex, Ey, Ez;
    pGrid Bx, By, Bz;

    schnek::MPICartSubdivision<Grid> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> E_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> B_parameters;
    schnek::pParametersGroup spaceVars;

    Vector initE;
    Vector initB;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void registerData();
    void fillValues();
  public:
    void init();
    void execute();

    static Index getGlobalMax() { return globalMax; }
    const schnek::DomainSubdivision<Grid> &getSubdivision() { return subdivision; };
};

#endif // MPULSE_MPULSE_H
