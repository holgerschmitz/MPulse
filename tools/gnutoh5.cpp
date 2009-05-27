#include <schnek/matrix.h>

#include "../src/hdfstream.h"
#include "../src/mpulse.h"

#include <sstream>
#include <fstream>

int main(int argc, char**argv) {
  char *infile = argv[1];
  char *outfile = argv[2];
  
  int X, Y, Z;
  int x, y, z;
  double val;
  
  DataGrid grid;

  std::ifstream input(infile);
  input >> X >> Y >> Z;
  grid.resize(GridIndex(X,Y,Z));
  
  for (int i=0; i<X; ++i)
    for (int j=0; j<Y; ++j)
      for (int k=0; k<Z; ++k)
      {
        input >> x >> y >> z >> val;
        grid(i,j,k) = val;
      }
  
  MatrixContainer<double, 3, schnek::MatrixAssertCheck> gridContainer;
  
  gridContainer.active = true;
  gridContainer.grid = &grid;
  gridContainer.global_min = GridIndex(0,0,0);
  gridContainer.global_max = GridIndex(X-1,Y-1,Z-1);
  
  
  HDFostream output(outfile);
  if (argc>4) {
    output.setBlockName(argv[3]);
  }
  output << gridContainer;
  output.close();


}
