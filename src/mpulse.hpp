/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_HPP_
#define MPULSE_HPP_


#include <schnek/grid.hpp>
#include <schnek/variables.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

#include <boost/ref.hpp>

#include <set>

#define GridChecker schnek::GridAssertCheck

typedef schnek::Array<int, 3> IndexType;
typedef schnek::Array<double, 3> Vector;

typedef schnek::Field<double, 3, GridChecker> Field;
typedef boost::shared_ptr<Field> pField;

typedef schnek::Range<int, 3> Range;
typedef schnek::Array<bool, 3> Stagger;

static const Stagger exStaggerYee(true,  false, false);
static const Stagger eyStaggerYee(false, true,  false);
static const Stagger ezStaggerYee(false, false, true);

static const Stagger bxStaggerYee(false, true,  true);
static const Stagger byStaggerYee(true,  false, true);
static const Stagger bzStaggerYee(true,  true,  false);

static const double PI     = 3.141592653589793238462643383279502884L;
static const double clight = 299792458;
static const double clight2 = clight*clight;
static const double mu_0 = 4e-7*PI;
static const double eps_0 = 1/(mu_0*clight2);


class FieldSolver;

class Simulation : public schnek::Block, public schnek::BlockContainer<FieldSolver> {
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
    schnek::Array<schnek::pParameter, 3> x_par, E_par, B_par;
    schnek::pParametersGroup spaceVars;
    schnek::Array<double, 3> initE, initB;

    double time;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
  public:
    Simulation();
    void init();
    void execute();

    static IndexType getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif /* MPULSE_HPP_ */
