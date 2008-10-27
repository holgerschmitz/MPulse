#include <schnek/matrix.h>

#include "../src/hdfstream.h"
#include "../src/mpulse.h"

#include <sstream>
#include <fstream>

class color
{
  public:
    color(unsigned char r_, unsigned char g_ , unsigned char b_)
      : r(r_), g(g_), b(b_) {}
    unsigned char r, g, b;
};

color getColor(double val)
{
  unsigned char r, g, b;
  
  if (val < 0) val = 0;
  
  double mval = 3*val;
  if (mval < 1)
  {
    return color(0,int(mval*255),255);
  }
  else if (mval < 2)
  {
    return color(255,255,int((2-mval)*255));
  }
  else if (mval < 3)
  {
    return color(int((mval-2)*255),int((3-mval)*255),0);
  }

  return color(255,0,0);

}

int main(int argc, char**argv) {
  char *coord = argv[1];
  std::string spos = argv[2];
  char *fname = argv[3];
  
  double amp = 1;
  
//  if (argc>4) {
//    std::string samp = argv[4];
//    std::istringstream ampstream(samp);
//    ampstream >> amp;
//  }
 
  DataGrid grid;
  HDFistream input(fname);
//  if (argc>5) {
//    input.setBlockName(argv[5]);
//  }
  input >> grid;
  input.close();
  
  int index[3];
  int first, second;
  std::istringstream posstream(spos);
  int pos;
  posstream >> pos;
  
  GridIndex low = grid.getLow();
  GridIndex high = grid.getHigh();
  
  switch (coord[0])
  {
    case 'x' : first = 1; second = 2; index[0] = pos; break;
    case 'y' : first = 0; second = 2; index[1] = pos; break;
    case 'z' : first = 0; second = 1; index[2] = pos; break;
  }
  
  int dim1 = high[first]-low[first] + 1;
  int dim2 = high[second]-low[second] + 1;
  std::cout << "P6\n" << dim1 << " " << dim2 <<"\n255\n";
  
  for (index[first]=low[first]; index[first]<=high[first]; ++index[first])
  {
    for (index[second]=low[second]; index[second]<=high[second]; ++index[second])
    {
      double val = grid(index[0], index[1], index[2]);
      color col = getColor(0.5*(amp*val+1));
      std::cout.put(col.r);
      std::cout.put(col.g);
      std::cout.put(col.b);
    }
  }
  
}
