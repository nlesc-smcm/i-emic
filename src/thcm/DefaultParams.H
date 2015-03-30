#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

//! in this class we provide static functions to set 
//! reasonable default parameters so that the user doesn't
//! need parameter lists with hundreds of entries. We choose
//! to do this in a central place. Anything that doesn't get
//! a default value here is left to the corresponding package/class,
//! e.g. LOCA, OceanPreconditioner etc. We override some of these defaults
//! because e.g. the preconditioner may choose default values that are
//! typically feasible on all machines, while for THCM it is advantageous
//! to use the best possible solvers by default (like MRILU).
class DefaultParams
{

public:

//! create the "THCM" sublist 
	static void THCM(Teuchos::ParameterList& list);

//! create the "NOX" sublist
	static void NOX(Teuchos::ParameterList& list);        

//! create the "LOCA" sublist
	static void LOCA(Teuchos::ParameterList& list);

//! create the "Time Stepping" sublist
	static void TimeStepping(Teuchos::ParameterList& list);
  
//! default settings for periodic orbit driver. The input list
//! should be the base list, changes are made to "THCM" and "NOX"
//! sublists.
	static void PeriodicOrbits(Teuchos::ParameterList& list);

};//class DefaultParams
