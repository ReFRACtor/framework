#include "stokes_coefficient_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(StokesCoefficient,
				 SubStateVectorArrayStokesCoefficient);

template<class Archive>
void StokesCoefficientImpBase::serialize(Archive & ar,
				const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayStokesCoefficient);
}

FP_IMPLEMENT(StokesCoefficientImpBase);
#endif
