#include "ils_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void IlsImpBase::serialize(Archive & ar,
			   const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Ils)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverSampleGrid)
    & FP_NVP_(sample_grid) & FP_NVP_(high_res_extension)
    & FP_NVP_(desc_band_name) & FP_NVP_(hdf_band_name);
}

template<class Archive>
void IdentityIls::serialize(Archive & ar,
			    const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IlsImpBase);
}

FP_IMPLEMENT(IlsImpBase);
FP_IMPLEMENT(IdentityIls);
#endif
