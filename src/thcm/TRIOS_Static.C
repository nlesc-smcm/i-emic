#define CLASSNAME "Static"
#include "TRIOS_Static.H"
#include "TRIOS_ftrace.H"
#include "TRIOS_Domain.H"
#include "Teuchos_RCP.hpp"
#include <mpi.h>

#ifdef DEBUGGING
#include "HYMLS_Tools.H"
#endif

namespace TRIOS
  {
  Teuchos::RCP<Domain> Static::domain_;
  
  void Trace(std::string s1, std::string s2, std::string s3)
    {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
#ifndef DEBUGGING    
    std::cerr << "PID "<< rank << " " << s1 << " " << s2 << "::"<<s3<<std::endl;
#else
  if (HYMLS::Tools::InitializedIO())
    {
    HYMLS::Tools::deb() << s1 << " " << s2 << "::" << s3 << std::endl; 
    }
  else
    {
    std::cerr << "PID "<< rank << " " << s1 << " " << s2 << "::"<<s3<<std::endl;
    }
#endif    
    }

  
  }

