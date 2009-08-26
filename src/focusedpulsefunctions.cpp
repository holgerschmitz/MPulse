#include "shortpulseinject.h"
#include <cmath>
#include "specfunc.h"

#define Power pow


double E = exp(1.);
double Pi = M_PI;

typedef FocusedPulseDataGenerator::Complex Complex;

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

inline Complex Sqrt(double x)
{
  return sqrt(Complex(x,0));
}

inline Complex Sqrt(Complex x)
{
  return sqrt(x);
}

FocusedPulseDataGenerator::Complex 
  FocusedPulseDataGenerator::FieldFuncs(double kx, double ky, double z, double t, int fieldid)
{

  double c = lightspeed;
  double p = Phase;
  
  double T = length;
  double T2 = T*T;
  double T4 = T2*T2;
  
  double Factor = (Power(Pi,1.5)*Power(Sigma,2)*\[Tau])/2.;
  double kp = sqrt(kx*kx+ky*ky);
  
  Complex Kph = -Sqrt(Power(k0,2) - Power(kp,2));

  Complex  Sph = -Power(sp,2) + s*spp)/Power(s,2);

  Complex  Tph = (-4*sp + Power(kp,2)*Power(s,2)*sp)/s + (Complex
(0,2)*(c*Sqrt (Power (Kph,2))*t + k0*z))/(c*Sqrt (Power (Kph,2))))/2.;

  Complex  phA = Sqrt (-(Power (kp,2)*Power (sp,2)) + 4*Sph - Power
(kp,2)*s*spp + (Complex (0,2)*Power (kp,2)*z)/(Power (c,2)*Power (Kph,3)) -
Power (Tau,2))/2.;

  Complex  phB = (Complex (0,0.5)*Tph)/phA;

  Complex  WA = Power (E,Power (Tau,2)*Power (Omega0,2))*Sqrt (Pi)*W
(-phB);

  Complex  WAM = Power (E,Power (Tau,2)*Power (Omega0,2))*Sqrt
(Pi)*W (phB);

  Complex  phC = 2*Dph*DTph - DDph*Tph;

  Complex  phD = Tph + Power (Tau,2)*Omega0;

  Complex  phE = -((phB*phD)/Tph);

  Complex  phF = -2*Dph*DTph + DDph*phD;

  Complex  phG = -Power(phD,2)/(4.*Dph) + Complex (0,2)*Phi;

  Complex  EphG = Power (E,phG);

  Complex  phH = (-Power(Tph,2) + 4*Dph*Power (Tau,2)*Power
(Omega0,2))/(4.*Dph);

  Complex  EphH = Power (E,phH);


  Complex  phI = EphG*phD + EphH*Tph;

  Complex  phJ = (-(Power (kp,2)*Power (s,2)) - Complex (0,4)*(Kph*z +
Phi) - 4*Power (Tau,2)*Power (Omega0,2))/4.;

  Complex  phK = 2*(EphG + EphH);

  Complex  phL = Dph*DTph - DDph*phD;

  Complex  phM = Dph*DTph - DDph*Tph;

  Complex  phN = phB/Tph;

  Complex  cfB1 = phK*Sqrt (Pi);

  Complex  cfB2 = Complex (0,1)*(Power (E,2*I*Phi) + Power (E,Power
(Tau,2)*Power (Omega0,2))) - 2*phI*phN*Sqrt (Pi);

  Complex  cfC2 = phN*Tph;

  Complex  cfD2 = phD*phN;
  

  Complex WA, WB;
  Complex SA, SB;

  if (imag(phB)>0)
  {
    SA = 2*Exp(Power(Tau,2)*Power(Omega0,2))*Sqrt(Pi)
    phJ -= phB*phB;
    WA = -Exp(phB*phB+Power(Tau,2)*Power(Omega0,2))*Sqrt(Pi)*Faddeeva_2(phB);
  }
  else
  {
    SA = 0.0;
    WA = Exp(Power(Tau,2)*Power(Omega0,2))*Sqrt(Pi)*Faddeeva_2(-phB);
  }
  
  if (imag(phE)<0)
  {
    SB = 2.0*Exp(2.0*I*Phi)*Sqrt(Pi)
    phJ -= phE*phE;
    WB = -Exp(phE*phE+2*I*Phi)*Sqrt(Pi)*Faddeeva_2(-phE);
  }
  else
  {
    SB = 0.0;
    WB = Exp(2*I*Phi)*Sqrt(Pi)*Faddeeva_2(phE);
  }
  
  Complex Psi = -I*Exp(phJ)*phN*(cfB1 - (WA+SA) - (WB+SB));
  Complex Psit = 2*Exp(phJ)*Power(phN,2)*(cfB2 + cfC2*(WA+SA) + cfD2*(WB+SB));
  Complex Psiz = (Exp(phJ)*Power(phN,2)*(cfB3 - cfC3*(WA+SA) + cfD3*(WB+SB)))/Dph;
  Complex Psizz =  -((Exp(phJ)*Power(phN,4)*(cfB5 - cfC5*(WA+SA) - cfD5*(WB+SB)))
    /Power(Dph,2));
  Complex Psitz = -(Exp(phJ)*phN*(cfB4 + cfC4*(WA+SA) + cfD4*(WB+SB)))
    /(2.*Power(Dph,2));
    
  Ex = -Factor*Psitz;
  Ey = -Factor*YComp*Psitz;
  Ez =  Factor*(kx+YComp*ky)*Psit;
  
  Bx = Factor*(-(kx+YComp*ky)*ky*Psi - YComp*Psizz);
  By = Factor*( (kx+YComp*ky)*kx*Psi + Psizz);
  Bz = Factor*  (ky+YComp*kx)*Psiz;
}


#undef Power
#undef Sqrt
