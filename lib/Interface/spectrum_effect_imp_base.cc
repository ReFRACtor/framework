#include "spectrum_effect_imp_base.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
SUB_STATE_VECTOR_ARRAY_SERIALIZE(SpectrumEffect, SubStateVectorArraySpectrumEffect);

template<class Archive>
void SpectrumEffectImpBase::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArraySpectrumEffect);
}

FP_IMPLEMENT(SpectrumEffectImpBase);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SpectrumEffectImpBase, SpectrumEffect)
REGISTER_LUA_END()
#endif
