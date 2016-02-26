#include "CouplingBlock.H"
#include "GlobalDefinitions.H"

//------------------------------------------------------------------
// constructor
CouplingBlock::CouplingBlock(std::string const &coupling)
	:
	coupling_(coupling)
{
	
}

//------------------------------------------------------------------
void CouplingBlock::initialize()
{
	
}

//------------------------------------------------------------------
void CouplingBlock::applyMatrix(SuperVector const &v, SuperVector &out)
{
	if (coupling_.compare("OA") == 0)
		applyOA(v, out);
	else if (coupling_.compare("AO") == 0)
		applyAO(v, out);
	else
	{
		WARNING("Invalid coupling mode!", __FILE__, __LINE__);
	}
}

//------------------------------------------------------------------
void CouplingBlock::applyOA(SuperVector const &v, SuperVector &out)
{
	
}

//------------------------------------------------------------------
void CouplingBlock::applyAO(SuperVector const &v, SuperVector &out)
{
}
