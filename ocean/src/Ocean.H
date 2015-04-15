#ifndef OCEAN_H
#define OCEAN_H

#include "GlobalDefinitions.H"
// For now this is a dummy class

class Ocean
{
public:
	//constructor
	Ocean()
		{
			//INFO("Ocean constructor called");
			std::cout << "Ocean constructor called" << std::endl;
		}
	//destructor
	~Ocean()
		{
			//INFO("Ocean destructor called");
		}
};
#endif
