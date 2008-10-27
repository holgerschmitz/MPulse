#include "pulseinit.h"
#include "globals.h"
#include "storage.h"



void GaussPulseInit::init(Storage &fields)
{
  std::cerr << "=== INITIALIZING Gaussian Pulse ===\n";
  double ex = by;
  double ey = -bx;
  
  double dx = Globals::instance().gridDX();
  double dy = Globals::instance().gridDY();
  double dz = Globals::instance().gridDZ();
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  
  GridIndex gridLow = Globals::instance().gridLow();
  GridIndex gridHigh = Globals::instance().gridHigh();

  DataGrid &Ex = fields.getGrid("Ex");
  DataGrid &Ey = fields.getGrid("Ey");
  DataGrid &Ez = fields.getGrid("Ez");
  DataGrid &Bx = fields.getGrid("Bx");
  DataGrid &By = fields.getGrid("By");
  DataGrid &Bz = fields.getGrid("Bz");

  GridIndex low = fields.getLow();
  GridIndex high = fields.getHigh();
  
  double nkz = dz*kz;
  double rr0 = r0*r0;
  double zz0 = z0*z0;

  double xh = 0.5*double(gridHigh[0] + gridLow[0]);
  double yh = 0.5*double(gridHigh[1] + gridLow[1]);

  
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
    {
//      double iie = i-xh;
//      double jje = j-yh;
//      double rre = (iie*iie*dx2 + jje*jje*dy2) / rr0;
      
//      double iib = i+0.5-xh;
//      double jjb = j+0.5-yh;
//      double rrb = (iib*iib*dx2 + jjb*jjb*dy2) / rr0;
      
      for (int k=low[2]+bound; k<=high[2]-bound; ++k)
      {
        double nze = k*dz - zc;
        double nzb = (k+0.5)*dz - zc;
        
        double kze = k*nkz;
        double kzb = (k+0.5)*nkz;
        
 //       double ampe = cos(kze - C*rre)*exp(-nze*nze/zz0 - rre);
 //       double ampb = cos(kzb - C*rrb)*exp(-nzb*nzb/zz0 - rrb);
        double ampe = cos(kze)*exp(-nze*nze/zz0);
        double ampb = cos(kzb)*exp(-nzb*nzb/zz0);
        
        Ex(i,j,k) = ex*ampe;
        Ey(i,j,k) = ey*ampe;
        Ez(i,j,k) = 0;
        
        Bx(i,j,k) = bx*ampb;
        By(i,j,k) = by*ampb;
        Bz(i,j,k) = 0;
      }
    }
  
}


ParameterMap* GaussPulseInit::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["k"] = WParameter(new ParameterValue<double>(&kz,1));
  (*pm)["r0"] = WParameter(new ParameterValue<double>(&r0,1));
  (*pm)["z0"] = WParameter(new ParameterValue<double>(&z0,1));
  (*pm)["zc"] = WParameter(new ParameterValue<double>(&zc,1));
  (*pm)["C"] = WParameter(new ParameterValue<double>(&C,1));
  (*pm)["Bx"] = WParameter(new ParameterValue<double>(&bx,1));
  (*pm)["By"] = WParameter(new ParameterValue<double>(&by,1));
  (*pm)["bound"] = WParameter(new ParameterValue<int>(&bound,0));
  
  return pm;
}
