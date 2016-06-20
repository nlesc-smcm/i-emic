//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPO_H
#define TOPO_H

#include "TopoDecl.H"
#include "GlobalDefinitions.H"

#include <vector>
#include <sstream>

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
Topo(Model model, ParameterList pars)
	:
	model_(model),
	pars_(pars),

	nMasks_(pars->get("Number of mask files", 0))	
{
	INFO("Topo constructor...");
	loadMaskFiles();
	loadMaskArrays();	
	INFO("Topo constructor... done");
}

//==================================================================
// constructor
template<typename Model, typename ParameterList>
Topo<Model, ParameterList>::
~Topo()
{
	INFO("Topo destructor");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMaskFiles()
{
	INFO(" Topo: loading mask files...");
	
	landMaskFiles_ = std::vector<std::string>(nMasks_);
	std::stringstream fileID;
	for (int i = 0; i != nMasks_; ++i)
	{
		fileID << "Mask file " << i;		
		landMaskFiles_[i] = pars_->get(fileID.str().c_str(), "");
		assert(landMaskFiles_[i] != "");
		INFO("  " << fileID.str() << ": " << landMaskFiles_[i]);
		fileID.str("");
		fileID.clear();
	}
	INFO(" Topo: loading mask files... done");
}

//==================================================================
template<typename Model, typename ParameterList>
void Topo<Model, ParameterList>::loadMaskArrays()
{
	INFO(" Topo: loading mask arrays...");
	landMaskArrays_ = std::vector<LandMask>(nMasks_);
	
	for (int i = 0; i != nMasks_; ++i)
		landMaskArrays_[i] = model_->getLandMask(landMaskFiles_[i]);
	INFO(" Topo: loading mask arrays... done");
}


#endif
