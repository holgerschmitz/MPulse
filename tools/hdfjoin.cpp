#include <hdf5.h>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <new>
#include <schnek/grid.hpp>


typedef schnek::Array<int, 3> Pos;

struct Bounds
{
  Bounds() {}
  Bounds(Pos low_, Pos high_) : low(low_), high(high_) {}
  Bounds(const Bounds &b) : low(b.low), high(b.high) {}
  Pos low;
  Pos high;
  int *Low() { return &(low[0]); }
  int *High() { return &(high[0]); }
};

struct hdffileinfo
{
  hid_t file_id;
  hid_t sid;
  hid_t dataset;
};

std::vector<Bounds> boundaries;
Bounds fullbound;

void init(int np)
{
  boundaries.resize(np);
  fullbound.low = Pos(1000,1000,1000);
  fullbound.high = Pos(0,0,0);
}

void readlimits(int numproc, std::string boundfile)
{
  // reformat boundfile string  
  size_t ppos = boundfile.find("%p");
  if (ppos != std::string::npos) boundfile.replace(ppos, 2, "%1$03d");
  
  for (int pr=0; pr<numproc; ++pr)
  {
    Pos low(0,0,0);
    Pos high(0,0,0);
    
    std::string fname = str(boost::format(boundfile) % pr);
    std::cerr << " | " << fname << "\n";
    
    std::ifstream boundf(fname.c_str());
    std::string line;
    
    while (std::getline(boundf, line).good())
    {
      std::istringstream linestr(line);
      std::string token;
      linestr >> token;
      if (token == "Low:")
        linestr >> low[0] >> low[1] >> low[2];
      if (token == "High:")
        linestr >> high[0] >> high[1] >> high[2];
      
    }
    
    for (int k=0; k<3; ++k)
    {
      if (low[k] <fullbound.low[k])  fullbound.low[k] =low[k];
      if (high[k]>fullbound.high[k]) fullbound.high[k]=high[k];
    }
    boundaries[pr] = Bounds(low,high);
  }
  
  std::cerr << " | Full extent: (" 
    << fullbound.low[0] <<","<<fullbound.low[1] <<","<< fullbound.low[2]<<") -- ("
    << fullbound.high[0] <<","<<fullbound.high[1] <<","<< fullbound.high[2]<<")\n";
}

hdffileinfo createhdf(std::string fname)
{
  hsize_t dims[3];
  for (int k=0; k<3; ++k) dims[k] = fullbound.high[k]-fullbound.low[k]+1;
  hdffileinfo info;
  
  info.file_id 
    = H5Fcreate (fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
  info.sid 
    = H5Screate_simple (3, dims, NULL);
    
  info.dataset 
    = H5Dcreate(info.file_id, "data", H5T_NATIVE_DOUBLE, info.sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  
  return info;
}

void writechunks(int numproc, hdffileinfo output, std::string infile)
{
  // reformat infile string  
  size_t ppos = infile.find("%p");
  if (ppos != std::string::npos) infile.replace(ppos, 2, "%1$03d");
  //std::cout << infile << std::endl;
  
  for (int pr=0; pr<numproc; ++pr)
  {
    std::string fname = str(boost::format(infile) % pr);
    
//    std::cout << fname << std::endl;
//    continue;
  
    std::cerr << " | " << fname << "\n";
    
    hsize_t dims[3];
    hsize_t start[3];
    for (int k=0; k<3; ++k)
    {
      dims[k] = boundaries[pr].high[k]-boundaries[pr].low[k]+1;
      start[k] = boundaries[pr].low[k];
    }
    
    std::cerr << "Dimensions:\n";
    std::cerr << "low:" << boundaries[pr].low[0] << " "
      << boundaries[pr].low[1] << " "<< boundaries[pr].low[2]  << "\n";
    std::cerr << "high:" << boundaries[pr].high[0] << " "
      << boundaries[pr].high[1] << " "<< boundaries[pr].high[2]  << "\n";
    
    double *data = new (std::nothrow) double[dims[0]*dims[1]*dims[2]];
    
    if (data==0)
    {
      std::cerr << "Could not allocate memory of size " << dims[0]*dims[1]*dims[2]/(1024*1024.) << "MB\n";
      exit(-1);
    }
    
    std::cerr << "Allocated memory of size " << 8*dims[0]*dims[1]*dims[2]/(1024*1024.) << "MB\n";
    for (int i=0; i<dims[0]*dims[1]*dims[2]; ++i)
    {
      data[i] = 0.001*i;
    }
    std::cerr << "Data pointer " << data << "\n";
    
    
    // Reading from input file
    std::cerr << " | | opening ... ";
    hid_t file_id = H5Fopen (fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset = H5Dopen(file_id, "data",H5P_DEFAULT);
    std::cerr << "reading ... ";
    H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    std::cerr << "closing ... ";
    H5Dclose(dataset);
    H5Fclose (file_id);
    std::cerr << "done\n";
    
    // writing chunk to output file
    std::cerr << " | | writing chunk\n";
    hid_t file_dataspace = H5Dget_space(output.dataset);
    
    H5Sselect_hyperslab(file_dataspace,  H5S_SELECT_SET, 
                        start, NULL, dims, NULL);
    hid_t mem_dataspace = H5Screate_simple(3, dims, NULL);
    H5Dwrite(output.dataset, 
             H5T_NATIVE_DOUBLE, 
             mem_dataspace, 
             file_dataspace,	    
             H5P_DEFAULT, 
             data);
    
    H5Sclose(mem_dataspace);
    H5Sclose(file_dataspace);

    std::cerr << " | | deleting data\n";
    delete[] data;
  }
  
}

void closehdf(hdffileinfo output)
{
  H5Dclose(output.dataset);
  H5Sclose(output.sid);
  H5Fclose(output.file_id);
}

namespace po = boost::program_options;

int main(int argc, char **argv)
{
  int numproc;
  std::string boundfile;
  std::string infile;
  std::string outfile;
  
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("np", po::value<int>(&numproc), 
          "set number of processes")
      ("input", po::value<std::string>(&infile), 
          "set input file name base (use %p for the processor number)")
      ("output", po::value<std::string>(&outfile), 
          "set output file name")
      ("bfile", po::value<std::string>(&boundfile)->default_value("boundary%p.dat"), 
          "set boundary file name");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if ((vm.count("help")>0) || 
      (vm.count("np")==0) || 
      (vm.count("output")==0)) {
      std::cerr << desc << "\n";
      return 1;
  }

  std::cerr << "Number of processes: " << numproc << std::endl;
  
  std::cerr << "Initializing\n";
  init(numproc);
  std::cerr << "Reading Limits\n";
  readlimits(numproc, boundfile);
  std::cerr << "Creating Output File\n";
  hdffileinfo globalfile = createhdf(outfile);
  std::cerr << "Writing Chunks\n";
  writechunks(numproc, globalfile, infile);
  std::cerr << "Closing Output File\n";
  closehdf(globalfile);
  std::cerr << "Done\n";
}
