#include <iostream>
#include <fstream>
#include <schnek/grid.hpp>

int sizex = 4514;
int sizey = 2030;

double maxval = 1.5;

typedef schnek::Grid<double, 2, schnek::GridAssertCheck> DataGrid;
typedef DataGrid::IndexType Index;

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

int main(int argc, char **argv)
{
    if (argc<3)
    {
	std::cerr << "Usage gnutopnm <infile> <outfile>\n";
	exit(-1);
    }

    std::ifstream input(argv[1]);

    DataGrid field(Index(0,0), Index(sizex, sizey));

    int x, y;
    double val;
    double max = 0.0, min=0.0;
    int cnt = 0;

    input >> x >> y >> val;

    
    while (! input.eof())
    { 
	if ((x>=0) && (y>=0) && (x<=sizex) && (y<=sizey)) 
	{
	  field(x,y) = val; 
	  if (val>max) max=val;
	  if (val<min) min=val;
	  ++cnt;
	}
	input >> x >> y >> val;
	//std::cout << x << y << val << std::endl;
    }

    std::cerr << "Read " << cnt << " values. Maximum = " << max  << " Minimum = " << min << std::endl;
    input.close();

    std::ofstream output(argv[2]);

    output << "P6\n" << sizey+1 << " " << sizex+1 <<"\n255\n";
  
    //double amp = 1/maxval;
    for (int i=sizex; i>=0; --i)
    {
	for (int j=0; j<=sizey; ++j)
	{
	    double val = field(i,j);
	    color col = getColor((val-min)/(max-min));
	    output.put(col.r);
	    output.put(col.g);
	    output.put(col.b);
	}
    }
    output.close();
}
