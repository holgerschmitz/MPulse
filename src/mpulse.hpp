/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/grid.hpp>
#include <schnek/variables/block.hpp>
#include <schnek/variables/blockcontainer.hpp>

#define GridChecker schnek::GridAssertCheck
#define ArrayChecker schnek::ArrayAssertArgCheck

typedef schnek::Array<int, 3, ArrayChecker> IndexType;
typedef schnek::Array<double, 3, ArrayChecker> Vector;
typedef schnek::Field<double, 3, GridChecker> Field;
typedef boost::shared_ptr<Field> pField;
typedef schnek::Range<int, 3, ArrayChecker> Range;
typedef schnek::Array<bool, 3, ArrayChecker> Stagger;

static const double clight = 299792458.0;
static const double clight2 = clight*clight;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true );

static const Stagger bxStaggerYee(false, true,  true );
static const Stagger byStaggerYee(true,  false, true );
static const Stagger bzStaggerYee(true,  true,  false);

class FDTDSolver;

class Simulation : public schnek::Block, public schnek::BlockContainer<FDTDSolver> {
  private:
    static Simulation *instance;
    IndexType globalMax;
    IndexType gridSize;
    Vector size;
    Vector dx;

    double cflFactor;
    double dt;

    double tMax;

    schnek::MPICartSubdivision<Field> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, 3> x_parameters;
    schnek::Array<schnek::pParameter, 3> E_parameters;
    schnek::Array<schnek::pParameter, 3> B_parameters;
    schnek::pParametersGroup spaceVars;

    Vector initE;
    Vector initB;

    Range innerRange;
    pField Ex, Ey, Ez;
    pField Bx, By, Bz;
    void fillValues();
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void registerData();
  public:
    Simulation();
    void init();
    void execute();

    static IndexType getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif // MPULSE_MPULSE_H
