#include <schnek/matrix.h>

#include "../src/hdfstream.h"
#include "../src/mpulse.h"

#include <sstream>
#include <fstream>

int main(int argc, char**argv) {
  char *fname = argv[1];
  
  DataGrid2d grid;
  HDFistream input(fname);
  
  input >> grid;
  input.close();
    
  GridIndex2d low = grid.getLow();
  GridIndex2d high = grid.getHigh();
  
  
  for (int i=low[0]; i<=high[0]; ++i)
  {
    for (int j=low[1]; j<=high[1]; ++j)
    {
      std::cout << i << " " << j << " "  << grid(i,j) << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
