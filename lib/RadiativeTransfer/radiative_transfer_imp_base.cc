#include "radiative_transfer_imp_base.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
SUB_STATE_VECTOR_ARRAY_SERIALIZE(RadiativeTransferRetrievable, SubStateVectorArrayRadiativeTransferRetrievable);

template<class Archive>
void RadiativeTransferImpBase::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayRadiativeTransferRetrievable);
}

FP_IMPLEMENT(RadiativeTransferImpBase);
#endif
