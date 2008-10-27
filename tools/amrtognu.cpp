#include <schnek/matrix.h>

#include "../src/hdfstream.h"

#include <sstream>
#include <fstream>

int main(int argc, char**argv) {
  char *fname = argv[1];
    
  DataGrid grid;
  HDFistream input(fname);
  input.setBlockName("data");
 
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
  
  for (index[first]=low[first]; index[first]<=high[first]; ++index[first])
  {
    for (index[second]=low[second]; index[second]<=high[second]; ++index[second])
    {
      std::cout << index[first] << " " << index[second] << " " 
        << grid(index[0], index[1], index[2]) << std::endl;
    }
    std::cout << std::endl;
  }
  
}
