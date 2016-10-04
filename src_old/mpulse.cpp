#include "globals.h"
#include "process.h"

#ifndef SINGLE_PROCESSOR
#define MPICH_SKIP_MPICXX
// extern "C" {
#include <mpi.h>
// }
#endif

#include <fstream>
#include <string>
#include <unistd.h>

// Print some initial diagnostics and start the simulation.
//---------------------------------------------------------------
//begin of main
int main (int argc, char** argv) {
  // if SINGLE_PROCESSOR is not defined, the MPI protocol is used 

#ifndef SINGLE_PROCESSOR
  MPI_Init(&argc, &argv);
#endif

//  std::cerr << "Waiting for debugger\n";
//  int DebugWait = 1;
//  while(DebugWait);

  // set  the static member of parameters
  Globals::setArgc(argc);
  // set  the static member of parameters
  Globals::setArgv(argv);

  //--------------------------------------------------------------- 
  // if specified get setup parameters from file
  // default setup file name, to be parsed if specified differently in the
  // argument
  std::string setupfilename = "setup.dat";
  if (argc > 1) setupfilename = argv[1];
  
  //buffer for current working directory
  char buf[200];
  //get current working directory
  getcwd(buf,200);
  std::cerr << "*******************************************************************************\n";
  std::cerr << "****************************** This is MPulse *********************************\n";
  std::cerr << "*******************************************************************************\n";
  std::cerr << "CWD is " << buf << "\n";
  
  //open setup file
  std::ifstream setupfile(setupfilename.c_str());
  std::cerr << "MAIN: Reading Inputfile: " << setupfilename.c_str() << " ...\n";
  //if "setup.dat" does not exist, display error message
  if (!setupfile)
  {
    std::cerr << "MAIN: Could not open: " << setupfilename.c_str() << " ...\n"; 
    exit(-1);
  }

  Process *process = new Process();
  process->Rebuild(setupfile);
  
  std::cerr << "                          *** MAIN:  INITIALIZING ***\n";
  
  //turn over control to process class, initialize...
  process->init();
  std::cout << "                             *** MAIN: RUNNING ***\n";
  
  //.. and run the process
  process->run();
  
  //delete the allocated process
  delete process;
  
  std::cerr << "*************************** MPulse has finished ***************************\n";
  return 0;
}

// end of main
