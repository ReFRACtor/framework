#include "aerosol_extinction_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(AerosolExtinction,
                                 SubStateVectorArrayAerosolExtinction);

template<class Archive>
void AerosolExtinctionImpBase::serialize(Archive & ar,
                                         const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayAerosolExtinction)
    & FP_NVP(press) 
    & FP_NVP_(aerosol_name);
}

FP_IMPLEMENT(AerosolExtinctionImpBase);
#endif

ArrayAd<double, 1> AerosolExtinctionImpBase::aerosol_extinction
(Pressure::PressureGridType Gtype) const
{
  fill_cache();
  if(Gtype == Pressure::NATIVE_ORDER ||
     (int) Gtype == (int) press->type_preference())
    return aext;
  return ArrayAd<double, 1>(aext.value().reverse(blitz::firstDim),
			    aext.jacobian().reverse(blitz::firstDim));
}

AutoDerivative<double>
AerosolExtinctionImpBase::extinction_for_layer(int i) const
{
  ArrayAd<double, 1> a = aerosol_extinction();
  range_check(i, 0, a.rows() - 1); 
  return (a(i) + a(i + 1)) / 2;
}
