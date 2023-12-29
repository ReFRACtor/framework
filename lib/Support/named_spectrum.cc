#include "named_spectrum.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NamedSpectrum::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Spectrum)
    & FP_NVP_(name) & FP_NVP_(index);
}

FP_IMPLEMENT(NamedSpectrum);
FP_OBSERVER_SERIALIZE(NamedSpectrum);
FP_OBSERVER_SERIALIZE(NamedSpectrumPtr);
FP_OBSERVER_SERIALIZE(NamedSpectrumPtrVec);
#endif

