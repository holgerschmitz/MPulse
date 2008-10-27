#include <schnek/matrix.h>

#include "../src/hdfstream.h"
#include "../src/mpulse.h"

#include <sstream>
#include <fstream>

int main(int argc, char**argv) {
  char *coord = argv[1];
  std::string spos = argv[2];
  char *fname = argv[3];
  
  DataGrid grid;
  HDFistream input(fname);
  if (argc>4) {
    input.setBlockName(argv[4]);
  }
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
    case 'x' : 
	first = 1; 
	second = 2; 
	index[0] = pos; 
	if ((pos<low[0]) || (pos>high[0])) return 0;
	break;
    case 'y' :
	first = 0; 
	second = 2; 
	index[1] = pos;
	if ((pos<low[1]) || (pos>high[1])) return 0;
	break;
    case 'z' : 
	first = 0;
	second = 1;
	index[2] = pos;
	if ((pos<low[2]) || (pos>high[2])) return 0;
	break;
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
  std::cout << std::endl;
}
