#include "observation_sv_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(ObservationSv, SubStateVectorArrayObservationSv);

template<class Archive>
void ObservationSvImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayObservationSv);
}

FP_IMPLEMENT(ObservationSvImpBase);
#endif
