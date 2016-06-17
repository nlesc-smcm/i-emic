//======================================================================
// Author T.E.Mulder --> t.e.mulder@uu.nl
//======================================================================
#ifndef TOPODECL_H
#define TOPODECL_H

//! This class manages a Model that allows changing its topography. A
//! continuation with Topo<Model> model is then possible through the
//! continuation parameter delta_.
//!
//! We expect a Model, which should be a (reference counted) pointer
//! to an I-EMIC model. A Model should maintain its own Vector,
//! which we expect to be of pointer type as well and named VectorPtr.
//!
//! The templated type ParameterList should be a pointer 
//! to a parameter storing object having a .get() method.

template<typename Model, typename ParameterList>
class Continuation
{
	using Vector = typename Model::element_type::VectorPtr;
	
	Model model_;
	
	
};


#endif
