#include "temperature_level.h"

#include "linear_interpolate.h"
#include "fp_serialize_support.h"

#include <boost/bind/bind.hpp>
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;
using namespace boost::placeholders;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void TemperatureLevel::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TemperatureImpBase);
}

FP_IMPLEMENT(TemperatureLevel);

#endif

//-----------------------------------------------------------------------
/// Set up Temperature
//-----------------------------------------------------------------------

TemperatureLevel::TemperatureLevel(const blitz::Array<double, 1> Temp,
                                   const boost::shared_ptr<Pressure>& Press,
                                   boost::shared_ptr<StateMapping> Map)
{
    init(Temp, Press, Map);
}

//-----------------------------------------------------------------------
/// This calculates temperature grid to use for layer retrieval.
//-----------------------------------------------------------------------

void TemperatureLevel::calc_temperature_grid() const
{
  blitz::Array<double, 1> press_profile( pressure_profile() );
  ArrayAd<double, 1> mapped_state_coeff = mapping->mapped_state(coeff);

  if (press_profile.rows() != mapped_state_coeff.rows()) {
    std::stringstream err_msg;
    err_msg << "Size of pressure grid: "
            << press_profile.rows()
            << " != size of temperature levels: "
            << mapped_state_coeff.rows();
    throw Exception(err_msg.str());
  }

  std::vector<AutoDerivative<double> > plist;
  std::vector<AutoDerivative<double> > tlist;
  for(int i = 0; i < press_profile.rows(); ++i) {
    tlist.push_back(mapped_state_coeff(i));
    plist.push_back(press_profile(i));
  }

  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > lin_type;
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
std::string TemperatureLevel::state_vector_name_i(int i) const
{
  return "Temperature (Kelvin) for Press Lvl " +
    boost::lexical_cast<std::string>(i + 1);
}

boost::shared_ptr<Temperature> TemperatureLevel::clone() const
{
  return boost::shared_ptr<TemperatureLevel>
    (new TemperatureLevel(coeff.value(), press->clone(), mapping->clone()));
}

void TemperatureLevel::print(std::ostream& Os) const
{
  Os << "Temperature Level\n";
}
