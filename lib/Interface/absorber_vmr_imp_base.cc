#include "absorber_vmr_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(AbsorberVmr, SubStateVectorArrayAbsorberVmr);

template<class Archive>
void AbsorberVmrImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayAbsorberVmr)
    & FP_NVP_(gas_name);
}

FP_IMPLEMENT(AbsorberVmrImpBase);
#endif
