/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include "types.hpp"

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>
#include <boost/enable_shared_from_this.hpp>


class FieldSolver;
class EMFields;

class MPulse : public schnek::Block,
               public schnek::BlockContainer<FieldSolver>,
               public schnek::BlockContainer<EMFields>,
               public boost::enable_shared_from_this<MPulse>
{
  private:
    static MPulse *instance;
    Index globalMax;
    Index gridSize;
    Vector size;
    Vector dx;

    double cflFactor;
    double dt;


    double tMax;
    double time;
    Range innerRange;
    pField pEx, pEy, pEz;
    pField pBx, pBy, pBz;

    schnek::MPICartSubdivision<Field> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> E_parameters;
    schnek::Array<schnek::pParameter, DIMENSION> B_parameters;
    schnek::pParametersGroup spaceVars;

    Vector initE;
    Vector initB;
  private:
    void initFields();
    void fillValues();
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
  public:
    MPulse();
    void init();
    void execute();

    static Index getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static double getDt() { return instance->dt; }
    static Vector getSize() { return instance->size; }
    static double getTime() { return instance->time; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif // MPULSE_MPULSE_H
