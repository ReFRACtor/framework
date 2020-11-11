#include "meteorology.h"
#include "fp_serialize_support.h"
#include "log_interpolate.h"
#include "old_constant.h"
#include <boost/algorithm/string.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Meteorology::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(Meteorology);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(Meteorology);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

typedef Array<double, 1> (Meteorology::*f_val)() const;
typedef Array<double, 1> (Meteorology::*f_interp)(const Array<double, 1>&) const;

typedef Array<double, 1> (Meteorology::*f_vmr)(const std::string&) const;
typedef Array<double, 1> (Meteorology::*f_vmr_interp)(const std::string&, const Array<double, 1>&) const;

REGISTER_LUA_CLASS(Meteorology)
.def("pressure_levels", &Meteorology::pressure_levels)
.def("specific_humidity", ((f_val) &Meteorology::specific_humidity))
.def("specific_humidity", ((f_interp) &Meteorology::specific_humidity))
.def("vmr", ((f_vmr) &Meteorology::vmr))
.def("vmr", ((f_vmr_interp) &Meteorology::vmr))
.def("temperature", ((f_val) &Meteorology::temperature))
.def("temperature", ((f_interp) &Meteorology::temperature))
.def("surface_pressure", &Meteorology::surface_pressure)
.def("windspeed", &Meteorology::windspeed)
.def("windspeed_u", &Meteorology::windspeed_u)
.def("windspeed_v", &Meteorology::windspeed_v)
REGISTER_LUA_END()
#endif

Array<double, 1> Meteorology::specific_humidity(const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(specific_humidity(), Pressure_level);
}

Array<double, 1> Meteorology::vmr(const std::string& Species, const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(vmr(Species), Pressure_level);
}

Array<double, 1> Meteorology::temperature(const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(temperature(), Pressure_level);
}

double Meteorology::windspeed() const
{
    return sqrt( sqr(windspeed_u()) + sqr(windspeed_v()) ); 
}

blitz::Array<double, 1> Meteorology::vmr(const std::string& Species) const
{
    std::string species_upper = Species;
    boost::to_upper(species_upper);
    if (species_upper == "H2O") {
        return h2o_vmr();
    } else {
        Exception err;
        err << "Can not return VMR species " << Species << " handling has not been defined.";
        throw err;
    }   
}

blitz::Array<double, 1> Meteorology::h2o_vmr() const
{
    Array<double, 1> s = specific_humidity();
    Array<double, 1> vmr(s.shape());
    vmr = s / (1 - s) * OldConstant::molar_weight_dry_air / OldConstant::molar_weight_water;
    return vmr;
}

Array<double, 1> Meteorology::interpolate_to_grid(const Array<double, 1>& Profile, const Array<double, 1>& Dest_pressure_levels) const
{
    LogLogInterpolate<double, double> prof_interp(pressure_levels().begin(), pressure_levels().end(), Profile.begin());
    Array<double, 1> prof_res(Dest_pressure_levels.shape());
    for(int i = 0; i < Dest_pressure_levels.rows(); ++i) {
        prof_res(i) = prof_interp(Dest_pressure_levels(i));
    }
    return prof_res;
}

