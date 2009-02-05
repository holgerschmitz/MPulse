#include "shortpulseinject.h"
#include <cmath>
#include "specfunc.h"

ShortPulseInjectSourceFunc::Complex 
  ShortPulseInjectSourceFunc::Exfunc(double x, double y, double z, double t)
{
  double c = lightspeed;
  double p = Phase;
  double tb = 0.5*om0*t;
  
  double T = length;
  double T2 = T*T;
  double T4 = T2*T2;
  
  double r2 = x*x+y*y;
  Complex I(0,1), I2(0,2), I4(0,4);
  Complex q(z,ZRl);
  Complex rT = ZRl/q;
  Complex rT4T = rT*T4;
  Complex rpiT = sqrt(M_PI)*rT;
  Complex rpT = -rT*rT/ZRl;

  Complex eph2T = exp(2.*I*(p + 2.*tb));

  Complex tpT = t - z/c - r2/(2.*c*q);
  Complex tpT2 = tpT*tpT;
  Complex tpT3 = tpT*tpT2;
  Complex tpT4 = tpT*tpT3;
  
  Complex tpzT = - 1/c - r2/(2.*c*q*q);
  
  Complex anT = - T*tb/tpT + I*tpT/T;
  Complex apT =   T*tb/tpT + I*tpT/T;
  
  Complex anT2 = anT*anT;
  Complex apT2 = apT*apT;
  
  Complex eanT2 = exp(anT2);
  Complex eapT2 = exp(apT2);

  Complex wp = Faddeeva_2(apT);
  Complex wn = Faddeeva_2(anT);
 
  return (exp(-I*p - I2*tb - 
       tpT2/T2)*(eanT2*rpiT*
        (rpT*T2*tpT3 + (-(apT2*rT4T) + rT*T2*tpT2 - 
             2.*rT*tpT4)*tpzT)*wn - 
       I*(anT*eapT2*eph2T*rpT*rT4T*tpT - 
          I2*rpiT*rpT*T2*tpT3 - 
          I2*eph2T*rpiT*rpT*T2*tpT3 + 
          anT*eanT2*rT*rT4T*tpzT + 
          apT2*eanT2*rT*rT4T*tpzT + 
          2.*anT*apT2*eanT2*rT*rT4T*tpzT + 
          anT2*eapT2*eph2T*rT*rT4T*tpzT - 
          I2*rpiT*rT*T2*tpT2*tpzT - 
          I2*eph2T*rpiT*rT*T2*tpT2*tpzT - 
          4.*anT*eapT2*eph2T*rT*rT*T2*tpT2*tpzT + 
          I4*rpiT*rT*tpT4*tpzT + 
          I4*eph2T*rpiT*rT*tpT4*tpzT + 
          apT*((1. + 2.*anT2)*eapT2*eph2T*rT*rT4T*tpzT + 
             eanT2*(rpT*rT4T*tpT - 
                4.*rT*rT*T2*tpT2*tpzT)) + 
          I*eapT2*eph2T*rpiT*
           (rpT*T2*tpT3 + 
             (-(anT2*rT4T) + rT*T2*tpT2 - 2.*rT*tpT4)*tpzT)*
           wp)))/(rpiT*T4*tpT2);
}

ShortPulseInjectSourceFunc::Complex
  ShortPulseInjectSourceFunc::Bfunc(double x, double y, double z, double t, bool bx)
{
  double c = lightspeed;
  double p = Phase;
  double tb = 0.5*om0*t;
  
  double T = length;
  double T2 = T*T;
  double T4 = T2*T2;
  
  double r2 = x*x+y*y;
  Complex I(0,1), I2(0,2), I4(0,4), nI(0,-1), nI2(0,-2);
  Complex q(z,ZRl);
  Complex rT = ZRl/q;
  Complex rT4T = rT*T4;
  Complex rpiT = sqrt(M_PI)*rT;
  Complex rpT = -rT*rT/ZRl;

  Complex eph2T = exp(2.*I*(p + 2.*tb));

  Complex tpT = t - z/c - r2/(2.*c*q);
  Complex tpT2 = tpT*tpT;
  Complex tpT3 = tpT*tpT2;
  Complex tpT4 = tpT*tpT3;
  
  Complex tpzT = - 1/c - r2/(2.*c*q*q);
  Complex tpzT2 = tpzT*tpzT;
  Complex tpxT = -x/(c*q);
  Complex tpyT = -y/(c*q);
  Complex tpxT2 = tpxT*tpxT;
  Complex tpyT2 = tpyT*tpyT;
  
  Complex anT = - T*tb/tpT + I*tpT/T;
  Complex apT =   T*tb/tpT + I*tpT/T;
  
  Complex anT2 = anT*anT;
  Complex apT2 = apT*apT;
  
  Complex eanT2 = exp(anT2);
  Complex eapT2 = exp(apT2);

  Complex wp = Faddeeva_2(apT);
  Complex wn = Faddeeva_2(anT);
 
  Complex BfA = (exp(nI*p - I2*tb - 
       tpT2/T2)*(rT*T2*tpT*tpxyT*
        (nI*apT*eanT2*rT*T2 - 
          I*anT*eapT2*eph2T*rT*T2 - 
          2*rpiT*tpT2 - 2*eph2T*rpiT*tpT2 + 
          eanT2*rpiT*tpT2*wn + 
          eapT2*eph2T*rpiT*tpT2*wp) + 
       tpxT*tpyT*
        (eanT2*rpiT*(-(apT2*rT4T) + rT*T2*tpT2 - 2*rT*tpT4)*
           wn - I*
           (rT*(anT*eanT2*rT4T + apT2*eanT2*rT4T + 
                2*anT*apT2*eanT2*rT4T + 
                anT2*eapT2*eph2T*rT4T + 
                apT*eapT2*eph2T*rT4T + 
                2*anT2*apT*eapT2*eph2T*rT4T - 
                I2*rpiT*T2*tpT2 - 
                I2*eph2T*rpiT*T2*tpT2 - 
                4*apT*eanT2*rT*T2*tpT2 - 
                4*anT*eapT2*eph2T*rT*T2*tpT2 + 
                I4*rpiT*tpT4 + 
                I4*eph2T*rpiT*tpT4) - 
             I*eapT2*eph2T*rpiT*
              (anT2*rT4T - rT*T2*tpT2 + 2*rT*tpT4)*wp)))
     )/(rpiT*T4*tpT2)
 
  Complex BfB = (exp(nI*p - I2*tb - 
       tpT2/T2)*(rT*T2*tpT*tpxxT*
        (nI*apT*eanT2*rT*T2 - 
          I*anT*eapT2*eph2T*rT*T2 - 
          2*rpiT*tpT2 - 2*eph2T*rpiT*tpT2 + 
          eanT2*rpiT*tpT2*wn + 
          eapT2*eph2T*rpiT*tpT2*wp) + 
       tpxT2*
        (eanT2*rpiT*(-(apT2*rT4T) + rT*T2*tpT2 - 2*rT*tpT4)*
           wn - I*
           (rT*(anT*eanT2*rT4T + apT2*eanT2*rT4T + 
                2*anT*apT2*eanT2*rT4T + 
                anT2*eapT2*eph2T*rT4T + 
                apT*eapT2*eph2T*rT4T + 
                2*anT2*apT*eapT2*eph2T*rT4T - 
                I2*rpiT*T2*tpT2 - 
                I2*eph2T*rpiT*T2*tpT2 - 
                4*apT*eanT2*rT*T2*tpT2 - 
                4*anT*eapT2*eph2T*rT*T2*tpT2 + 
                I4*rpiT*tpT4 + 
                I4*eph2T*rpiT*tpT4) - 
             I*eapT2*eph2T*rpiT*
              (anT2*rT4T - rT*T2*tpT2 + 2*rT*tpT4)*wp)))
     )/(rpiT*T4*tpT2)
 
  Complex BfC = (exp(nI*p - I2*tb - 
       tpT2/T2)*(rT*T2*tpT*tpyyT*
        (nI*apT*eanT2*rT*T2 - 
          I*anT*eapT2*eph2T*rT*T2 - 
          2*rpiT*tpT2 - 2*eph2T*rpiT*tpT2 + 
          eanT2*rpiT*tpT2*wn + 
          eapT2*eph2T*rpiT*tpT2*wp) + 
       tpyT2*
        (eanT2*rpiT*(-(apT2*rT4T) + rT*T2*tpT2 - 2*rT*tpT4)*
           wn - I*
           (rT*(anT*eanT2*rT4T + apT2*eanT2*rT4T + 
                2*anT*apT2*eanT2*rT4T + 
                anT2*eapT2*eph2T*rT4T + 
                apT*eapT2*eph2T*rT4T + 
                2*anT2*apT*eapT2*eph2T*rT4T - 
                I2*rpiT*T2*tpT2 - 
                I2*eph2T*rpiT*T2*tpT2 - 
                4*apT*eanT2*rT*T2*tpT2 - 
                4*anT*eapT2*eph2T*rT*T2*tpT2 + 
                I4*rpiT*tpT4 + 
                I4*eph2T*rpiT*tpT4) - 
             I*eapT2*eph2T*rpiT*
              (anT2*rT4T - rT*T2*tpT2 + 2*rT*tpT4)*wp)))
     )/(rpiT*T4*tpT2)
 
  Complex BfD = (exp(nI*p - I2*tb - 
       tpT2/T2)*(nI2*
        (I*(1 + eph2T)*rpiT*tpT2*
           (rppT*rT4T - 2*rT*
              (2*rpT*T2*tpT*tpzT + 
                rT*(T2*tpzT2 - 
                   2*tpT2*tpzT2 + T2*tpT*tpzzT))) + 
          rT*rT*T2*
           (apT2*eanT2*rT*T2*tpzT2 + 
             anT2*eapT2*eph2T*rT*T2*tpzT2 + 
             apT*((1 + 2*anT2)*eapT2*eph2T*rT*T2*
                 tpzT2 + 
                eanT2*tpT*
                 (2*rpT*T2*tpzT - 4*rT*tpT*tpzT2 + 
                   rT*T2*tpzzT)) + 
             anT*((1 + 2*apT2)*eanT2*rT*T2*tpzT2 + 
                eapT2*eph2T*
                 (2*rpT*T2*tpT*tpzT - 
                   4*rT*tpT2*tpzT2 + rT*T2*tpT*tpzzT
                   )))) - 
       eanT2*rpiT*(rppT*rT4T*tpT2 - 
          2*rT*(2*rpT*T2*tpT3*tpzT - 
             apT2*rT4T*tpzT2 + 
             rT*T2*tpT2*tpzT2 - 
             2*rT*tpT4*tpzT2 + rT*T2*tpT3*tpzzT))*
        wn - eapT2*eph2T*rpiT*
        (rppT*rT4T*tpT2 - 
          2*rT*(2*rpT*T2*tpT3*tpzT - 
             anT2*rT4T*tpzT2 + 
             rT*T2*tpT2*tpzT2 - 
             2*rT*tpT4*tpzT2 + rT*T2*tpT3*tpzzT))*
        wp))/
   (2.*rpiT*rT*T4*tpT2);
   
   if (bx) return BfA+BfC-YComp*BfD;
   else return BfD-BfB-BfA;
}
