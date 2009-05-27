#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>

#define Power pow
#define Sqrt sqrt

double Pi = M_PI;
#ifndef M_2PI
#define M_2PI       6.2831853071795864769252867665590   // 2Pi
#endif

std::ofstream debug1("debug1.out");
std::ofstream debug2("debug2.out");
std::ofstream debug3("debug3.out");
std::ofstream debug4("debug4.out");
std::ofstream debug5("debug5.out");
std::ofstream debug6("debug6.out");

typedef std::complex<double> Complex;

class ShortPulseInjectSourceFunc
{
  public:

    void setParam(double length_, 
                  double width_,
                  double om0_,
                  double Shift_,
                  double Phase_,
                  double amp_, 
                  double eps_);

    Complex Efunc(double x, double y, double z, double t);
    Complex Bfunc(double x, double y, double z, double t, bool bx);

  private:
    
    double lightspeed;

    Complex YComp;

    double length; // corresponds to pulse duration
    double width;
    double om0;
    double Shift;
    double Phase;

    double ZRl;
    
    double amp;
    double eps;
    
};

Complex Faddeeva_2(const Complex& z)
{
    // table 1: coefficients for h = 0.5  
    static double n1[12] =
             { 0.25, 1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0,
               20.25, 25.0, 30.25, 36.0 };
    static double e1[12] =
             { 0.7788007830714049,    0.3678794411714423,
               1.053992245618643e-1,  1.831563888873418e-2,
               1.930454136227709e-3,  1.234098040866795e-4,
               4.785117392129009e-6,  1.125351747192591e-7,
               1.605228055185612e-9,  1.388794386496402e-11,
               7.287724095819692e-14, 2.319522830243569e-16 };

    // table 2: coefficients for h = 0.53 
    static double n2[12] =
             { 0.2809, 1.1236, 2.5281, 4.4944, 7.0225, 10.1124,
               13.7641, 17.9776, 22.7529, 28.09, 33.9889, 40.4496 };
    static double e2[12] =
             { 0.7551038420890235,    0.3251072991205958, 
               7.981051630007964e-2,  1.117138143353082e-2,
               0.891593719995219e-3,  4.057331392320188e-5,
               1.052755021528803e-6,  1.557498087816203e-8,
               1.313835773243312e-10, 6.319285885175346e-13,
               1.733038792213266e-15, 2.709954036083074e-18 };
    
   // tables for Pade approximation 
   static double C[7] =
            { 65536.0, -2885792.0, 69973904.0, -791494704.0,
              8962513560.0, -32794651890.0, 175685635125.0 };
   static double D[7] =
            { 192192.0, 8648640.0, 183783600.0, 2329725600.0,
              18332414100.0, 84329104860.0, 175685635125.0 };


    double *n,*e,t,u,r,s,d,f,g,h;
    Complex c,d2,v,w,zz;
    int i;
    
    // use Pade approximation 
    s = norm(z);
    if (s < 1e-7) {
        zz = z*z;
        v  = exp(zz);
        c  = C[0];
        d2 = D[0];
        for (i = 1; i <= 6; i++) {
            c  = c  * zz + C[i];
            d2 = d2 * zz + D[i];
        }
        w = 1.0 / v + Complex(0.0,M_2_SQRTPI) * c/d2 * z * v;
        return w;

    // use trapezoid rule 
    } else {

        // select default table 1 
        n = n1;
        e = e1;
        r = M_1_PI * 0.5;
 
        // if z is too close to a pole select table 2 
        if (fabs(imag(z)) < 0.01 && fabs(real(z)) < 6.01) {
            h = modf(2*fabs(real(z)),&g);
            if (h < 0.02 || h > 0.98) {
                n = n2;
                e = e2;
                r = M_1_PI * 0.53;
            }
           }
        
        d = (imag(z) - real(z)) * (imag(z) + real(z));
        f = 4 * real(z) * real(z) * imag(z) * imag(z);

        g = h = 0.0;
        for (i = 0; i < 12; i++) {
            t = d + n[i];
            u = e[i] / (t * t + f);
            g += (s + n[i]) * u;
            h += (s - n[i]) * u;
        }
        u = 1 / s;
        c = r * Complex(imag(z) * (u + 2.0 * g),
                                real(z) * (u + 2.0 * h) );
        
        if (imag(z) < M_2PI) {
            s = 2.0 / r;
            t = s * real(z);
            u = s * imag(z);
            s = sin(t);
            h = cos(t);
            f = exp(- u) - h;
            g = 2.0 * exp(d-u) / (s * s + f * f);
            u = 2.0 * real(z) * imag(z);
            h = cos(u);
            t = sin(u);
            c += g * Complex( (h * f - t * s), -(h * s + t * f));
        }
        return c;
    }
}


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


void ShortPulseInjectSourceFunc::setParam(double length_,
                                          double width_,
                                          double om0_,
                                          double Shift_,
                                          double Phase_,
                                          double amp_, 
                                          double eps_)
{
  std::cerr << "ShortPulseInjectSourceFunc::setParam\n";
  length = length_;
  width  = width_;
  om0     = om0_;
  Shift  = Shift_;
  Phase  = 2*M_PI*Phase_;
    
  amp = amp_;
  eps = eps_;
  
  lightspeed = sqrt(1/eps);
  ZRl = 0.5*om0*width*width/lightspeed;

  YComp = Complex(0.0,0.0);
}

Complex ShortPulseInjectSourceFunc::Efunc(double x, double y, double z, double t)
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
           
  return -I*Ef/c;
}

Complex ShortPulseInjectSourceFunc::Bfunc(double x, double y, double z, double t, bool bx)
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

int main()
{
  ShortPulseInjectSourceFunc source;
  source.setParam(1e6, 2, 2, 0, 0, 1, 1);
  
  double y=0;
  double z=0;
  for (double x=-5; x<=5.01; x+=0.05)
  {
    for (double t=0; t<=8; t+=0.025)
    {
      Complex Ex = source.Efunc(x,y,z,t);
      Complex Bx = source.Bfunc(x,y,z,t,true);
      Complex By = source.Bfunc(x,y,z,t,false);
      debug1 << x << " " << t << " " << Ex.real() << "\n";
      debug2 << x << " " << t << " " << Ex.imag() << "\n";
      debug3 << x << " " << t << " " << Bx.real() << "\n";
      debug4 << x << " " << t << " " << Bx.imag() << "\n";
      debug5 << x << " " << t << " " << By.real() << "\n";
      debug6 << x << " " << t << " " << By.imag() << "\n";
    }
    debug1 << "\n";
    debug2 << "\n";
    debug3 << "\n";
    debug4 << "\n";
    debug5 << "\n";
    debug6 << "\n";
  }
}
