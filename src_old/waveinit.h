#ifndef MPULSE_WAVEINIT_H
#define MPULSE_WAVEINIT_H

#include "init.h"

class PlaneWaveInit : public FieldSimInit
{
  public:
    void init(Storage &fields);
  protected:
    /// build parametermap
    ParameterMap* MakeParamMap (ParameterMap* pm = NULL);
  private:
    int kx, ky, kz;
    double bx, by, bz;
    double eps;
};


#endif
