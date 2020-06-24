#include <boost/bind.hpp>
#include "absorber_vmr_level.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberVmrLevel::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbsorberVmrImpBase);
}

FP_IMPLEMENT(AbsorberVmrLevel);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevel, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<double, 1>&,
			  const blitz::Array<bool, 1>&,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrLevel::AbsorberVmrLevel
(const boost::shared_ptr<Pressure>& Press,
 const blitz::Array<double, 1>& Vmr, 
 const blitz::Array<bool, 1>& Vmr_flag,
 const std::string& Gas_name,
 boost::shared_ptr<Mapping> in_map)
{
  bool Mark_according_to_press = false;
  int Pdep_start = 0;
  init(Gas_name, Vmr, Vmr_flag, Press, Mark_according_to_press, Pdep_start, in_map);
}

AbsorberVmrLevel::AbsorberVmrLevel(const boost::shared_ptr<Pressure>& Press,
                                   const blitz::Array<double, 1>& Vmr,
                                   const bool Vmr_flag,
                                   const std::string& Gas_name,
                                   boost::shared_ptr<Mapping> in_map)
{
  bool Mark_according_to_press = false;
  int Pdep_start = 0;
  blitz::Array<bool, 1> flag(1);
  flag(0) = Vmr_flag;
  init(Gas_name, Vmr, flag, Press, Mark_according_to_press, Pdep_start, in_map);
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevel::clone() const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrLevel(press->clone(), coeff.value(),used_flag,
			  gas_name(), mapping->clone()));
}

blitz::Array<double, 1> AbsorberVmrLevel::pressure_profile() const
{
  return press->pressure_grid().value.value();
}

void AbsorberVmrLevel::calc_vmr() const
{
  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > vmrlist;
  ArrayAd<double, 1> fm_view_coeff = mapping->fm_view(coeff);
  for(int i = 0; i < press->pressure_grid().rows(); ++i) {
    vmrlist.push_back(fm_view_coeff(i));
    plist.push_back(press->pressure_grid()(i).value);
  }
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), vmrlist.begin()));
  vmr = boost::bind(&lin_type::operator(), lin, _1);
}

void AbsorberVmrLevel::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrLevel" + mapping->name() + ":\n"
     << "  Gas name:       " << gas_name() << "\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval Flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
