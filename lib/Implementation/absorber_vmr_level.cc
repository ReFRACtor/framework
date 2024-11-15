#include <boost/bind/bind.hpp>

#include "absorber_vmr_level.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;
using namespace boost::placeholders;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberVmrLevel::serialize(Archive & ar,
				 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbsorberVmrImpBase)
    & FP_NVP(coeff_pressure);
}

FP_IMPLEMENT(AbsorberVmrLevel);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
#include "state_mapping_at_indexes.h"

boost::shared_ptr<AbsorberVmr> absorber_vmr_level_create
(const boost::shared_ptr<Pressure>& Mapped_Press,
 const blitz::Array<double, 1>& Vmr,
 const blitz::Array<bool, 1>& Flag,
 const std::string& Gas_name)
{
  boost::shared_ptr<StateMapping> mapping =
    boost::make_shared<StateMappingAtIndexes>(Flag);

  boost::shared_ptr<AbsorberVmrLevel> abs_vmr = 
    boost::make_shared<AbsorberVmrLevel>(Mapped_Press, Vmr, Gas_name, mapping);

  return abs_vmr;
}

REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevel, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
     const blitz::Array<double, 1>&,
     const std::string&>())
.scope
[
 luabind::def("create", &absorber_vmr_level_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrLevel::AbsorberVmrLevel
(const boost::shared_ptr<Pressure>& Mapped_Press,
 const blitz::Array<double, 1>& Vmr,
 const std::string& Gas_name,
 boost::shared_ptr<StateMapping> in_map,
 const boost::shared_ptr<Pressure>& Coeff_Press)
{
  init(Gas_name, Vmr, Mapped_Press, in_map);

  if (!Coeff_Press) {
    coeff_pressure = Mapped_Press;
  } else {
    coeff_pressure = Coeff_Press;
  }

  if (coeff_pressure->pressure_grid().rows() != Vmr.rows()) {
    Exception err;
    err << "Coefficient pressure grid size: "
	<< coeff_pressure->pressure_grid().rows()
	<< " does not match VMR size: " << Vmr.rows()
	<< " for gas: " << Gas_name;
    throw err;
  }
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevel::clone() const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrLevel(mapped_pressure->clone(),
			  mapping->mapped_state(coeff.value()).value(),
			  gas_name(), mapping->clone(),
			  coeff_pressure->clone()));
}

blitz::Array<double, 1> AbsorberVmrLevel::pressure_profile() const
{
  return mapped_pressure->pressure_grid(Pressure::NATIVE_ORDER).value.value();
}

blitz::Array<double, 1> AbsorberVmrLevel::vmr_profile() const
{
  return mapping->mapped_state(coeff).value();
}

blitz::Array<double, 1> AbsorberVmrLevel::coeff_pressure_profile() const
{
  return coeff_pressure->pressure_grid(Pressure::NATIVE_ORDER).value.value();
}

void AbsorberVmrLevel::calc_vmr() const
{
  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > vmrlist;
  ArrayAd<double, 1> mapped_state_coeff = mapping->mapped_state(coeff);

  if(mapped_pressure->pressure_grid().rows() != mapped_state_coeff.rows()) {
    Exception err;
    err << "Mapped pressure grid size: " << mapped_pressure->pressure_grid().rows()
	<< " does not match VMR forward model grid size: " << mapped_state_coeff.rows();
    throw err;
  }

  for(int i = 0; i < mapped_pressure->pressure_grid().rows(); ++i) {
    vmrlist.push_back(mapped_state_coeff(i));
    plist.push_back(mapped_pressure->pressure_grid(Pressure::NATIVE_ORDER)(i).value);
  }

  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  if(mapped_pressure->type_preference() ==
     Pressure::PREFER_INCREASING_PRESSURE) {
    boost::shared_ptr<lin_type> lin
      (new lin_type(plist.begin(), plist.end(), vmrlist.begin()));
    vmr = boost::bind(&lin_type::operator(), lin, _1);
  } else {
    boost::shared_ptr<lin_type> lin
      (new lin_type(plist.rbegin(), plist.rend(), vmrlist.rbegin()));
    vmr = boost::bind(&lin_type::operator(), lin, _1);
  }
}


std::string AbsorberVmrLevel::state_vector_name_i(int coeff_idx) const
{
  blitz::Array<double, 1> press_grid(coeff_pressure->pressure_grid(Pressure::NATIVE_ORDER).value.value());
  std::stringstream sv_name;
  if(coeff.rows() == press_grid.rows()) {
    // Output the pressure associated with the retrieval value in millibars
    double press_val = coeff_pressure->pressure_grid(Pressure::NATIVE_ORDER)(coeff_idx).value.value();
    sv_name << gas_name() << " " << mapping->name() << " VMR at "
	    << std::fixed << std::setprecision(3) << (press_val / 100) << " hPa";
  } else {
    // If coefficient doesn't match the number of grid levels, we might
    // have something like StateMappingBasisMatrix. In any case, the
    // state vector isn't in terms of our pgrid, so use a more generic
    // description
    sv_name << gas_name() << " " << mapping->name() << " VMR at Level " << coeff_idx;
  }
  return sv_name.str();
}

void AbsorberVmrLevel::print(std::ostream& Os) const
{
  blitz::Array<double, 1> press_grid(coeff_pressure->pressure_grid(Pressure::NATIVE_ORDER).value.value());
  blitz::Array<double, 1> vmr_grid(vmr_profile());

  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrLevel\n"
     << "  Gas name: " << gas_name() << "\n"
     << "  StateMapping:  " << mapping->name() << "\n\n"
     << "      Pressure          VMR\n"
     << "  ------------ ------------\n";
  for(int i = 0; i < press_grid.rows(); i++) {
    Os << "  "
       << std::fixed << std::setprecision(3) << setw(8)
       << (press_grid(i) / 100) << " hPa" << " "
       << std::scientific << std::setprecision(5) << setw(12)
       << vmr_grid(i) << "\n";
  }
  opad.strict_sync();
}
