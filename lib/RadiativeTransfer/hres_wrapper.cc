#include "hres_wrapper.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void HresWrapper::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransferSingleWn)
    & FP_NVP_(rt);
}

FP_IMPLEMENT(HresWrapper);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<RadiativeTransfer>
hres_wrapper_create(const boost::shared_ptr<RadiativeTransfer>& Rt)
{
  boost::shared_ptr<RadiativeTransferSingleWn> rt =
    boost::dynamic_pointer_cast<RadiativeTransferSingleWn>(Rt);
  return boost::shared_ptr<RadiativeTransfer>
  (new HresWrapper(rt));
}
REGISTER_LUA_DERIVED_CLASS(HresWrapper, RadiativeTransfer)
.scope
[
 luabind::def("create", &hres_wrapper_create)
]
REGISTER_LUA_END()
#endif

// See base class for description
blitz::Array<double, 1> HresWrapper::stokes_single_wn
(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  return rt_->stokes_single_wn(Wn, Spec_index, Opt_prop);
}

// See base class for description
ArrayAd<double, 1> HresWrapper::stokes_and_jacobian_single_wn
(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop) const
{
  return rt_->stokes_and_jacobian_single_wn(Wn, Spec_index, Opt_prop);
}
