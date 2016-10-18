#include "shortpulseinject.hpp"
#include "specfunc.hpp"

#include <cmath>

#define Power pow
#define Sqrt sqrt

double E = exp(1.);
double Pi = M_PI;

typedef ShortPulseInjectSourceFunc::Complex Complex;

inline Complex operator*(int a, Complex b)
{
  return double(a)*b;
}

inline Complex operator*(Complex a, int b)
{
  return a*double(b);
}

inline Complex operator+(int a, Complex b)
{
  return double(a)+b;
}

ShortPulseInjectSourceFunc::Complex 
  ShortPulseInjectSourceFunc::Efunc(double x, double y, double z, double t)
{

  double c = lightspeed;
  double p = Phase;
  
  double T = length;
  double T2 = T*T;
  double T4 = T2*T2;
  double tb = 0.5*om0*T;
  
  double r2 = x*x+y*y;
  Complex I(0,1), I2(0,2), I4(0,4), nI(0,-1);
  Complex q(z,ZRl);
  Complex rT = ZRl/q;
  Complex rT4T = rT*T4;
  Complex rpiT = sqrt(M_PI)*rT;
  Complex rpT = -ZRl/(q*q);
  
  Complex eiphi = exp(2*I*p);

  Complex eph2T = exp(2.*I*(p + 2.*tb));

  Complex tpT = t - z/c - r2/(2.*c*q);
  Complex tpT2 = tpT*tpT;
  Complex tpT3 = tpT*tpT2;
  Complex tpT4 = tpT*tpT3;
  
  Complex tpzT = - 1/c - r2/(2.*c*q*q);
  
  Complex an = nI*tb + tpT/T;
  Complex ap =  I*tb + tpT/T;
 
  Complex A = exp(nI*p)/(rpiT*T4);
  Complex B = rT*T*(I*(1 + eiphi)*rpT*T2 - 
          rT*(T*(nI - 2*tb + 
                eiphi*(nI + 2*tb)) + 
             I2*(1 + eiphi)*tpT)*tpzT);
  Complex C = I*eiphi*rpiT*
        (T*tb + I*tpT)*
        (rpT*T2 + I2*rT*T*tb*tpzT - 
          2*rT*tpT*tpzT);
  Complex D = rpiT*(T*tb - I*tpT)*
        (I*rpT*T2 + 
          2*rT*(T*tb - I*tpT)*tpzT);
          
  Complex S = 0;
  
  if (imag(an)<0)
  {
    S += 2*A*C*exp(2*I*tpT*tb/T - tpT2/T2);
    an = -an;
    C = -C;
  }
  
  if (imag(ap)<0)
  {
    S += 2*A*D*exp(-2*I*tpT*tb/T - tpT2/T2);
    ap = -ap;
    D = -D;
  }
  
  Complex Ef = S + A*exp(-tb*tb)*( B + C*Faddeeva_2(an) - D*Faddeeva_2(ap));
           
  return -I*Ef;
}

ShortPulseInjectSourceFunc::Complex
  ShortPulseInjectSourceFunc::Bfunc(double x, double y, double z, double t, bool bx)
{
  double c = lightspeed;
  double p = Phase;
  
  double T = length;
  double T2 = T*T;
  double T4 = T2*T2;
  double tb = 0.5*om0*T;
  
  double r2 = x*x+y*y;
  Complex I(0,1), I2(0,2), I4(0,4), nI(0,-1), nI2(0,-2);
  Complex q(z,ZRl);
  Complex rT = ZRl/q;
  Complex rT2 = rT*rT;
  Complex rT4T = rT*T4;
  Complex rpiT = sqrt(M_PI)*rT;
  Complex rpT = -ZRl/(q*q);
  Complex rppT = 2.*ZRl/(q*q*q);

  Complex eiphi = exp(2*I*p);
  Complex eph2T = exp(2.*I*(p + 2.*tb));

  Complex tpT = t - z/c - r2/(2.*c*q);
  Complex tpT2 = tpT*tpT;
  Complex tpT3 = tpT*tpT2;
  Complex tpT4 = tpT*tpT3;
  
  Complex tpzT = - 1/c + r2/(2.*c*q*q);
  Complex tpzT2 = tpzT*tpzT;
  Complex tpxT = -x/(c*q);
  Complex tpyT = -y/(c*q);
  Complex tpxyT = 0.0;
  Complex tpxxT = -1./(c*q);
  Complex tpyyT = -1./(c*q);
  Complex tpzzT = -r2/(c*q*q*q);
  Complex tpxT2 = tpxT*tpxT;
  Complex tpyT2 = tpyT*tpyT;
  
  Complex aip = T*tb + I*tpT;
  Complex ain = T*tb - I*tpT;
  Complex aip2 = aip*aip;
  Complex ain2 = ain*ain;

  Complex an = nI*tb + tpT/T;
  Complex ap =  I*tb + tpT/T;
  
  Complex BfAA = exp(nI*p)/(Sqrt(Pi)*T4);
  Complex BfAB = -tpxT*tpyT*rT*T*(T*(nI - 2*tb + eiphi*(nI + 2*tb)) + I2*(1 + eiphi)*tpT)
        +I*T2*tpxyT*(1 + eiphi)*rT*T;
  Complex BfAC = -tpxT*tpyT*2*eiphi*rpiT*aip2
        +I*T2*tpxyT*eiphi*rpiT*(aip);
  Complex BfAD =-tpxT*tpyT*2*rpiT*(ain2)
        +I*T2*tpxyT*rpiT*(-(T*tb) + I*tpT);
        
  Complex BfAS = 0;
  
  Complex BfBA = exp(nI*p)/(Sqrt(Pi)*(T4));
  Complex BfBB =  - tpxT2*(rT*T*(T*(nI - 2*tb + eiphi*(nI + 2*tb)) + I2*(1 + eiphi)*tpT))
          + I*T2*tpxxT*(1 + eiphi)*rT*T;       
  Complex BfBC = - tpxT2*2*eiphi*rpiT*aip2
          + I*T2*tpxxT*eiphi*rpiT*(aip);     
  Complex BfBD = - tpxT2*2*rpiT*(ain2)
          + I*T2*tpxxT*rpiT*(-(T*tb) + I*tpT);
           
  Complex BfBS = 0;
 
  Complex BfCA = exp(nI*p)/(Sqrt(Pi)*(T4));
  Complex BfCB = - tpyT2*(rT*T*(T*(nI - 2*tb + eiphi*(nI + 2*tb)) + I2*(1 + eiphi)*tpT)) 
            + I*T2*tpxxT*(1 + eiphi)*rT*T;
  Complex BfCC = - tpyT2*2*eiphi*rpiT*aip2
            + I*T2*tpxxT*eiphi*rpiT*(aip);
  Complex BfCD = - tpyT2*2*rpiT*(ain2)
            + I*T2*tpxxT*rpiT*(-(T*tb) + I*tpT);
  
  Complex BfCS = 0;
  
  Complex BfDA =  exp(nI*p)/(2.*Sqrt(Pi)*T4);
  Complex BfDB =  2*T*(I2*(1 + eiphi)*rpT*T2*tpzT - 
            rT*(T*(nI - 2*tb + eiphi*(nI + 2*tb))*(tpzT2) + 
            I*(1 + eiphi)*(2*tpT*(tpzT2) - T2*tpzzT)));
  Complex BfDC = +(eiphi*rpiT*(rppT*rT4T + I2*rT*(aip)*
             (2*rpT*T2*tpzT + rT*(I2*T*tb*(tpzT2) - 2*tpT*(tpzT2) + T2*tpzzT))))/(rT2);
  Complex BfDD = +(rpiT*(rppT*rT4T - 2*rT*(ain)*(I2*rpT*T2*tpzT + 
              rT*(2*T*tb*(tpzT2) - I2*tpT*(tpzT2) + I*T2*tpzzT))))/(rT2);

  Complex BfDS = 0;

  if (imag(an)<0)
  {
    BfAS += 2*BfAA*BfAC*exp(2*I*tpT*tb/T - tpT2/T2);
    BfBS += 2*BfBA*BfBC*exp(2*I*tpT*tb/T - tpT2/T2);
    BfCS += 2*BfCA*BfCC*exp(2*I*tpT*tb/T - tpT2/T2);
    BfDS += 2*BfDA*BfDC*exp(2*I*tpT*tb/T - tpT2/T2);

    BfAC = -BfAC;
    BfBC = -BfBC;
    BfCC = -BfCC;
    BfDC = -BfDC;

    an = -an;
  }
  
  if (imag(ap)<0)
  {
    BfAS += 2*BfAA*BfAD*exp(-2*I*tpT*tb/T - tpT2/T2);
    BfBS += 2*BfBA*BfBD*exp(-2*I*tpT*tb/T - tpT2/T2);
    BfCS += 2*BfCA*BfCD*exp(-2*I*tpT*tb/T - tpT2/T2);
    BfDS += 2*BfDA*BfDD*exp(-2*I*tpT*tb/T - tpT2/T2);

    BfAD = -BfAD;
    BfBD = -BfBD;
    BfCD = -BfCD;
    BfDD = -BfDD;

    ap = -ap;
  }
        
  Complex wn = Faddeeva_2(an);
  Complex wp = Faddeeva_2(ap);

  Complex BfA = BfAS + BfAA*exp(-tb*tb)*(BfAB + BfAC*wn - BfAD*wp);
  Complex BfB = BfBS + BfBA*exp(-tb*tb)*(BfBB + BfBC*wn - BfBD*wp);
  Complex BfC = BfCS + BfCA*exp(-tb*tb)*(BfCB + BfCC*wn - BfCD*wp);
  Complex BfD = BfDS + BfDA*exp(-tb*tb)*(BfDB + BfDC*wn - BfDD*wp);

  
   if (bx) return I*(BfA + YComp*(BfC - BfD));
   else return I*(BfD - BfB - YComp*BfA);
}


#undef Power
#undef Sqrt
