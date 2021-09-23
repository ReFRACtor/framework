#include "absorber_vmr.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberVmr::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableAbsorberVmr);
}

FP_IMPLEMENT(AbsorberVmr);
FP_OBSERVER_SERIALIZE(AbsorberVmr);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AbsorberVmr)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<AbsorberVmr> >::*pbt1)(
        const std::vector<boost::shared_ptr<AbsorberVmr> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AbsorberVmr> >, VectorAbsorberVmr)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<AbsorberVmr> >::push_back))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return the vmr on the pressure grid.
//-----------------------------------------------------------------------

ArrayAd<double, 1> AbsorberVmr::vmr_grid(const Pressure& P,
					 Pressure::PressureGridType Gtype) const
{
  ArrayAd<double, 1> pgrid = P.pressure_grid(Gtype).convert(Unit("Pa")).value;
  blitz::Array<AutoDerivative<double>, 1> res(pgrid.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = volume_mixing_ratio(pgrid(i));
  return ArrayAd<double, 1>(res);
}
