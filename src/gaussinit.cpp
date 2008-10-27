#include "gaussinit.h"
#include "globals.h"
#include "storage.h"



void PlaneGaussPulseInit::init(Storage &fields)
{
  std::cerr << "=== INITIALIZING Gaussian Pulse ===\n";
  double ex = ky*bz-kz*by;
  double ey = kz*bx-kx*bz;
  double ez = kx*by-ky*bx;
  double bmag = sqrt(bx*bx + by*by + bz*bz);
  double factor = -bmag/sqrt(ex*ex + ey*ey + ez*ez);
  
  ex *= factor;
  ey *= factor;
  ez *= factor;
  
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
  
  double mx = 0.5*double(gridHigh[0] + gridLow[0]);
  double my = 0.5*double(gridHigh[1] + gridLow[1]);
  double mz = 0.5*double(gridHigh[2] + gridLow[2]);

  double nkx = kx/double(gridHigh[0] - gridLow[0] - 1);
  double nky = ky/double(gridHigh[1] - gridLow[1] - 1);
  double nkz = kz/double(gridHigh[2] - gridLow[2] - 1);  
  
  int offsetx = 0;
  int offsety = 0;
  int offsetz = 12;
  
  for (int i=low[0]+offsetx; i<=high[0]-offsetx; ++i)
    for (int j=low[1]+offsety; j<=high[1]-offsety; ++j)
      for (int k=low[2]+offsetz; k<=high[2]-offsetz; ++k)
      {
        double delta = 0.5;
        
/*
        double pose = (i-mx)*(i-mx)*nkx*nkx
                    + (j-my)*(j-my)*nky*nky 
                    + (k-mz)*(k-mz)*nkz*nkz;
        double posb = (i+delta-mx)*(i+delta-mx)*nkx*nkx 
                    + (j+delta-my)*(j+delta-my)*nky*nky 
                    + (k+delta-mz)*(k+delta-mz)*nkz*nkz;
        
        double ampe = exp(-pose);
        double ampb = exp(-posb);
*/
        double pose = (i-mx)*nkx
                    + (j-my)*nky 
                    + (k-mz)*nkz;
        double posb = (i+delta-mx)*nkx 
                    + (j+delta-my)*nky 
                    + (k+delta-mz)*nkz;
        
        double ampe = exp(-pose*pose);
        double ampb = exp(-posb*posb);
        
        Ex(i,j,k) = ex*ampe;
        Ey(i,j,k) = ey*ampe;
        Ez(i,j,k) = ez*ampe;
        
        Bx(i,j,k) = -bx*ampb;
        By(i,j,k) = -by*ampb;
        Bz(i,j,k) = -bz*ampb;
      }
  
}


ParameterMap* PlaneGaussPulseInit::MakeParamMap (ParameterMap* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["kx"] = WParameter(new ParameterValue<int>(&kx,1));
  (*pm)["ky"] = WParameter(new ParameterValue<int>(&ky,1));
  (*pm)["kz"] = WParameter(new ParameterValue<int>(&kz,1));
  (*pm)["Bx"] = WParameter(new ParameterValue<double>(&bx,1));
  (*pm)["By"] = WParameter(new ParameterValue<double>(&by,1));
  (*pm)["Bz"] = WParameter(new ParameterValue<double>(&bz,1));
  
  return pm;
}
