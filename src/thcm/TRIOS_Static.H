#ifndef _TRIOS_STATIC
#define _TRIOS_STATIC

#include "TRIOS_Domain.H"
#include "Teuchos_RCP.hpp"

namespace TRIOS 
{
  
	//! static functionality that makes programming some things easier.
	class Static 
    {
    public:
    
		//! register your domain with TRIOS (this allows constructing
		//! the block preconditioner and other classes in an Ifpack style,
		//! without passing in a domain object).
		static void SetDomain(Teuchos::RCP<Domain> domain)
			{
				domain_=domain;
			}
  
		//! if a domain has been registered using SetStaticDomain(), this returns
		//! the domain object. Otherwise, it returns Teuchos::null.
		static Teuchos::RCP<Domain> GetDomain()
			{
				return domain_;
			}
    
    private:
    
		//! domain object
		static Teuchos::RCP<Domain> domain_;
    };
}
#endif 
