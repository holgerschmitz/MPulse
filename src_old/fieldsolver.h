#ifndef MPULSE_FIELDSOLVER_H
#define MPULSE_FIELDSOLVER_H

#include "rebuild.h"
#include <cstdlib>

class Storage;
class Current;

class FieldSolver : public Rebuildable
{
  public:    
    virtual void initStorage(Storage *storage_) = 0;
    
    virtual void addCurrent(Current *current) 
    {
      std::cerr << "========================================================================\n";
      std::cerr << "=============================    ERROR     =============================\n";
      std::cerr << "========================================================================\n";
      std::cerr << "    Additional currents cannot be used when using this Field Solver     \n;";
      std::cerr << "========================================================================\n\n";
      exit(-1);
    }

    virtual void addMagCurrent(Current *current) 
    {
      std::cerr << "========================================================================\n";
      std::cerr << "=============================    ERROR     =============================\n";
      std::cerr << "========================================================================\n";
      std::cerr << "Additional magnetic currents cannot be used when using this Field Solver\n;";
      std::cerr << "========================================================================\n\n";
      exit(-1);
    }
    
    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;    
};

#endif
