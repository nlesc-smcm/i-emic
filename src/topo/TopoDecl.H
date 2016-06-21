//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPODECL_H
#define TOPODECL_H

#include <vector>

//! This class manages a Model that allows changing its topography. A
//! continuation with Topo<Model> model is then possible through the
//! continuation parameter delta_.
//!
//! We expect a Model, which should be a (reference counted) pointer
//! to an I-EMIC model. A Model should maintain its own Vector,
//! which we expect to be of pointer type as well and named VectorPtr.
//!
//! In addition to the requirements from the Continuation class, we
//! need a Model to implement member functions:
//!
//!   void          setLandMask()
//!   LandMask      getLandMask()
//! 
//! The templated type ParameterList should be a pointer 
//! to a parameter storing object having a .get() method.

template<typename Model, typename ParameterList>
class Topo
{
	using VectorPtr = typename Model::element_type::VectorPtr;
	using LandMask  = typename Model::element_type::LandMask;
	
	//! Model
	Model model_;

	//! Parameterlist
	ParameterList pars_;

	//! Number of masks we need to load
	int nMasks_;

	//! Vector containing all the filenames
	std::vector<std::string> landMaskFileNames_;

	//! Vector containing RCP's to the landmask arrays
	std::vector<LandMask> landMasks_;

	//! Continuation parameter
	double delta_;

public:
	//! default constructor
	Topo(){};
	//! constructor
	Topo(Model model, ParameterList pars);
	
	//! destructor
	~Topo();

private:
	//! load the mask filenames
	void loadMaskFileNames();

	//! load the masks
	void loadMasks();
};


#endif
