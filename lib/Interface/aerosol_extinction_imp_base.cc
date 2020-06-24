#include "aerosol_extinction_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(AerosolExtinction,
				 SubStateVectorArrayAerosolExtinction);

template<class Archive>
void AerosolExtinctionImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayAerosolExtinction)
    & FP_NVP_(aerosol_name);
}

FP_IMPLEMENT(AerosolExtinctionImpBase);
#endif
