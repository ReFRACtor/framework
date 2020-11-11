#include "aerosol_property_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(AerosolProperty,
				 SubStateVectorArrayAerosolProperty);

template<class Archive>
void AerosolPropertyImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayAerosolProperty);
}

FP_IMPLEMENT(AerosolPropertyImpBase);
#endif
