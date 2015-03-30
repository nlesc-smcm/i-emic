/**********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/

#ifndef SINGLETON_H
#define SINGLETON_H

#include "globdefs.H"
#include "Teuchos_RCP.hpp"

#if defined(GNU)||defined(G95)||defined(HUYGENS)
#define _INSTANTIATE_(type) template <typename type> type* Singleton<type>::instance;
#else
#define _INSTANTIATE_(type) type* Singleton<type>::instance=NULL;
#endif

template<typename T>
class Singleton
{
public:
	
	static T& Instance() {return *instance;}
  
	static bool Instantiated() {return instance != NULL;}
  
protected:
  
	Singleton(T* inst)
		{
			DEBUG("Singleton constructor...");
			if (instance != NULL)
			{
				INFO("Warning: re-instantiating singleton object");
				delete [] instance;
			}
			instance = inst;
		}
	
	virtual ~Singleton()
		{
			DEBUG("Destroy Singleton");
		}
	
	static T* instance;
	
};

#endif

