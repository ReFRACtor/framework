#include "radiative_transfer_retrievable.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void RadiativeTransferRetrievable::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransfer)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableRadiativeTransferRetrievable);
}

FP_IMPLEMENT(RadiativeTransferRetrievable);
FP_OBSERVER_SERIALIZE(RadiativeTransferRetrievable);
#endif
