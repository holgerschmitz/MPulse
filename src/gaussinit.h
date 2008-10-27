#ifndef MPULSE_GAUSSINIT_H
#define MPULSE_GAUSSINIT_H

#include "init.h"

class PlaneGaussPulseInit : public FieldSimInit
{
  public:
    void init(Storage &fields);
  protected:
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    int kx, ky, kz;
    double bx, by, bz;
};

#endif
