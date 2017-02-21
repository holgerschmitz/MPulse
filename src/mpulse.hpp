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

using namespace schnek;

static const double PI     = 3.141592653589793238462643383279502884L;
static const double clight = 299792458;
static const double clight2 = clight*clight;
static const double mu_0 = 4e-7*PI;
static const double eps_0 = 1/(mu_0*clight2);

static const int ghostCells = 2;

typedef Array<int, 2> IndexType;

extern Array<int, 2> globalMax;
extern MPICartSubdivision<Field<double, 2> > *subdivision;
extern Array<double, 2> dx;

class FieldSolver;

class SimulationBlock : public Block, public BlockContainer<FieldSolver> {
  private:
    MPICartSubdivision<Field<double, 2> > subdivision;
    IndexType gridSize;
    Array<double, 2> size;

    double cflFactor;
    double dt;
    double tMax;
    double time;
  protected:
    void initParameters(BlockParameters &parameters);
  public:
    void preInit();
    void execute();
};

#endif /* MPULSE_HPP_ */
