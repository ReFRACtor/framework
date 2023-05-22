#include "measured_radiance_field_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(MeasuredRadianceField, SubStateVectorArrayMeasuredRadianceField);

template<class Archive>
void MeasuredRadianceFieldImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayMeasuredRadianceField);
}

FP_IMPLEMENT(MeasuredRadianceFieldImpBase);
#endif
