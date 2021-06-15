#include "spectrum_sampling.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectrumSampling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpectrumSampling);
  ar & FP_NVP(nspectrometer);
}

template<class Archive>
void IdentitySpectrumSampling::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectrumSampling);
}

FP_IMPLEMENT(SpectrumSampling);
FP_IMPLEMENT(IdentitySpectrumSampling);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SpectrumSampling)
REGISTER_LUA_END()
#endif
