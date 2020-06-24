#include <boost/bind.hpp>
#include "absorber_vmr_fixed_level.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberVmrFixedLevel::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbsorberVmrImpBase)
    & FP_NVP(press_level);
}

FP_IMPLEMENT(AbsorberVmrFixedLevel);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrFixedLevel, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const boost::shared_ptr<PressureLevelInput>&,
			  const blitz::Array<bool, 1>&, 
			  const blitz::Array<double, 1>&,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrFixedLevel::AbsorberVmrFixedLevel
(const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<PressureLevelInput>& Press_level,	   
 const blitz::Array<bool, 1>& Used_flag, 
 const blitz::Array<double, 1>& Vmr,
 const std::string& Gas_name)
: AbsorberVmrImpBase(Gas_name, Vmr, Used_flag, Press),
  press_level(Press_level)
{
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrFixedLevel::clone() const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrFixedLevel(press->clone(), press_level, used_flag, coeff.value(),
			       gas_name()));
}

void AbsorberVmrFixedLevel::calc_vmr() const
{
  AutoDerivative<double> p = log(press->surface_pressure().value);
  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > vmrlist;
  for(int i = 0; i < press->pressure_grid().rows() - 1; ++i) {
    vmrlist.push_back(coeff(i));
    plist.push_back(press->pressure_grid()(i).value);
  }
  int i = press->pressure_grid().rows() - 1;
  double p1 = log(press_level->pressure_level()(i - 1));
  double p2 = log(press_level->pressure_level()(i));
  AutoDerivative<double> v = coeff(i - 1) + 
    (p - p1) * (coeff(i) - coeff(i - 1)) / (p2 - p1);
  plist.push_back(press->surface_pressure().value);
  vmrlist.push_back(v);
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> lin
    (new lin_type(plist.begin(), plist.end(), vmrlist.begin()));
  vmr = boost::bind(&lin_type::operator(), lin, _1);
}

void AbsorberVmrFixedLevel::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrFixedLevel:\n"
     << "  Gas name: " << gas_name() << "\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval Flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
