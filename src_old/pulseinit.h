#ifndef MPULSE_PULSEINIT_H
#define MPULSE_PULSEINIT_H

#include "init.h"

class GaussPulseInit : public FieldSimInit
{
  public:
    void init(Storage &fields);
  protected:
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    double kz;
    double r0;
    double z0;
    double zc;
    double C;
    double bx, by;
    int bound;
};

#endif
