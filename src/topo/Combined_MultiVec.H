#ifndef COMBINED_MULTIVEC
#define COMBINED_MULTIVEC

// This class combines two multivectors with 

class Combined_MultiVec
{
};

//------------------------------------------------------------------
// Specialization of MultiVectorTraits for Belos
//------------------------------------------------------------------
namespace Belos
{
	template<>
	class MultiVecTraits<double, Combined_MultiVec>
	{
		
	};

} // end of Belos namespace





#endif
