#include <boost/bind/bind.hpp>
#include "temperature_offset.h"
#include "linear_interpolate.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;
using namespace boost::placeholders;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void TemperatureOffset::serialize(Archive& ar,
                                  const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TemperatureImpBase);
}

FP_IMPLEMENT(TemperatureOffset);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureOffset, Temperature)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature Offset. 
//-----------------------------------------------------------------------

TemperatureOffset::TemperatureOffset(const boost::shared_ptr<Pressure>& Press,
                                     double Temp_offset)
{
  Array<double, 1> val(1);
  val(0) = Temp_offset;
  init(val, Press);
}

//-----------------------------------------------------------------------
/// This calculates temperature grid to use for layer retrieval. 
//-----------------------------------------------------------------------

void TemperatureOffset::calc_temperature_grid() const
{
  blitz::Array<double, 1> temp_profile( temperature_profile() );
  blitz::Array<double, 1> press_profile( pressure_profile() );

  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > tlist;
  if (press_profile.rows() != temp_profile.rows()) {
    std::stringstream err_msg;
    err_msg << "Size of pressure grid: "
            << press_profile.rows()
            << " != size of temperature levels: "
            << temp_profile.rows();
    throw Exception(err_msg.str());
  }
  for(int i = 0; i < press_profile.rows(); ++i) {
    AutoDerivative<double> t2 = temp_profile(i) + coefficient()(0);
    tlist.push_back(t2);
    plist.push_back(press_profile(i));
  }
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  if(plist.size() < 2)
    throw Exception("Must have at least 2 pressure level");
  if(plist[1].value() > plist[0].value()) {
    boost::shared_ptr<lin_type> lin(new lin_type(plist.begin(), plist.end(), tlist.begin()));
    cache.tgrid = boost::bind(&lin_type::operator(), lin, _1);
  } else {
    boost::shared_ptr<lin_type> lin(new lin_type(plist.rbegin(), plist.rend(), tlist.rbegin()));
    cache.tgrid = boost::bind(&lin_type::operator(), lin, _1);
  }
}

// See base class for description of this function.
std::string TemperatureOffset::state_vector_name_i(int UNUSED(i)) const
{
  return "Temperature Offset (Kelvin)";
}

//-----------------------------------------------------------------------
/// Uncertainty of temperature offset.
//-----------------------------------------------------------------------

double TemperatureOffset::temperature_offset_uncertainty() const
{
  if(cov.rows() == 0 ||
     cov(0,0) < 0)
    return 0.0;
  return sqrt(cov(0,0));
}

void TemperatureOffset::print(std::ostream& Os) const 
{
  Os << "TemperatureOffset\n";
}
