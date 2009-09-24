#include "waveinit.h"
#include "globals.h"
#include "storage.h"



void PlaneWaveInit::init(Storage &fields)
{
  std::cerr << "=== INITIALIZING Plane Wave ===\n";
  double ex = ky*bz-kz*by;
  double ey = kz*bx-kx*bz;
  double ez = kx*by-ky*bx;
  double bmag = sqrt(bx*bx + by*by + bz*bz);
  double factor = -bmag/sqrt(ex*ex + ey*ey + ez*ez);
  
  ex *= factor/sqrt(eps);
  ey *= factor/sqrt(eps);
  ez *= factor/sqrt(eps);
  
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
  
  double nkx = 2*M_PI*kx/double(gridHigh[0] - gridLow[0] - 1);
  double nky = 2*M_PI*ky/double(gridHigh[1] - gridLow[1] - 1);
  double nkz = 2*M_PI*kz/double(gridHigh[2] - gridLow[2] - 1);
  
  for (int i=low[0]; i<=high[0]; ++i)
    for (int j=low[1]; j<=high[1]; ++j)
      for (int k=low[2]; k<=high[2]; ++k)
      {
        int ii=i-1, jj=j-1, kk=k-1;
        double ampex = cos((ii+0.5)*nkx + jj*nky + kk*nkz);
        double ampey = cos(ii*nkx + (jj+0.5)*nky + kk*nkz);
        double ampez = cos(ii*nkx + jj*nky + (kk+0.5)*nkz);
        double ampbx = cos(ii*nkx + (jj+0.5)*nky + (kk+0.5)*nkz);
        double ampby = cos((ii+0.5)*nkx + jj*nky + (kk+0.5)*nkz);
        double ampbz = cos((ii+0.5)*nkx + (jj+0.5)*nky + kk*nkz);
        
        Ex(i,j,k) += ex*ampex;
        Ey(i,j,k) += ey*ampey;
        Ez(i,j,k) += ez*ampez;
        
        Bx(i,j,k) += bx*ampbx;
        By(i,j,k) += by*ampby;
        Bz(i,j,k) += bz*ampbz;
      }
  
}


ParameterMap* PlaneWaveInit::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["kx"] = WParameter(new ParameterValue<int>(&kx,1));
  (*pm)["ky"] = WParameter(new ParameterValue<int>(&ky,1));
  (*pm)["kz"] = WParameter(new ParameterValue<int>(&kz,1));
  (*pm)["Bx"] = WParameter(new ParameterValue<double>(&bx,1));
  (*pm)["By"] = WParameter(new ParameterValue<double>(&by,1));
  (*pm)["Bz"] = WParameter(new ParameterValue<double>(&bz,1));
  (*pm)["eps"] = WParameter(new ParameterValue<double>(&eps,1));
  
  return pm;
}
