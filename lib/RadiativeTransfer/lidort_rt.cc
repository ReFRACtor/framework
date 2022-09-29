#include "lidort_rt.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void LidortRt::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRt);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void LidortRt::save(Archive &ar,
		   const unsigned int UNUSED(version)) const
{
  // Save a little extra information if we aren't already saving full
  // state.
  if(!SpurrRt::serialize_full_state) {
    int number_stream_ = number_stream();
    int number_moment_ = number_moment();
    bool pure_nadir_ = rt_driver()->pure_nadir();
    bool do_multi_scatt_only_ = rt_driver()->do_multi_scatt_only();
    bool do_thermal_scattering_ = rt_driver()->do_thermal_scattering();
    ar & FP_NVP_(number_stream) & FP_NVP_(number_moment)
      & FP_NVP_(pure_nadir) & FP_NVP_(do_multi_scatt_only)
      & FP_NVP_(do_thermal_scattering);
  }
}

template<class Archive>
void LidortRt::load(Archive &ar,
		   const unsigned int UNUSED(version)) 
{
  if(!rt_driver_) {
    int number_stream_;
    int number_moment_;
    bool pure_nadir_, do_multi_scatt_only_, do_thermal_scattering_;
    ar & FP_NVP_(number_stream) & FP_NVP_(number_moment)
      & FP_NVP_(pure_nadir) & FP_NVP_(do_multi_scatt_only)
      & FP_NVP_(do_thermal_scattering);
    rt_driver_.reset(new LidortRtDriver
       (number_stream_, number_moment_, do_multi_scatt_only_, surface_type(),
	zen, pure_nadir_, do_solar_sources, do_thermal_emission,
	do_thermal_scattering_));
  }
}

FP_IMPLEMENT(LidortRt);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(LidortRt, RadiativeTransfer)
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&, 
                          bool, int, int, bool>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm The atmpsphere to use
/// \param Stokes_coef The Stokes coefficients.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir Flag for controlling azimuth dependence in the
///      output LIDORT will complain if user zenith is 0 and this is
///      not set, when not using ss correction mode
/// \param Number_streams Number of streams to use
/// \param Number_moments Number of moments to use in phase function (-1 to use
///     full resolution
/// \param Do_multi_scatt_only Whether to do multiple scattering calculations
///     only
//-----------------------------------------------------------------------

LidortRt::LidortRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                   const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                   const blitz::Array<double, 1>& Sza,
                   const blitz::Array<double, 1>& Zen,
                   const blitz::Array<double, 1>& Azm,
                   bool Pure_nadir,
                   int Number_streams,
                   int Number_moments,
                   bool Do_multi_scatt_only,
                   bool do_solar_sources,
                   bool do_thermal_emission,
                   bool do_thermal_scattering)
: SpurrRt(Atm, Stokes_coef, Sza, Zen, Azm, do_solar_sources, do_thermal_emission)
{   
  rt_driver_.reset(new LidortRtDriver(Number_streams, Number_moments, Do_multi_scatt_only, surface_type(), Zen, Pure_nadir, do_solar_sources, do_thermal_emission, do_thermal_scattering));
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void LidortRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "LidortRt\n";
  OstreamPad opad1(Os, "  ");
  SpurrRt::print(opad1, Short_form);
  opad1.strict_sync();
  OstreamPad opad(Os, "    ");
  Os << "  Number stream: " << number_stream() << "\n"
     << "  Number moments: " << number_moment() << "\n"
     << "  Use first order scattering calc: " 
     << (rt_driver()->do_multi_scatt_only() ? "True\n" : "False\n")
     << (rt_driver()->pure_nadir() ? "True\n" : "False\n");

  opad << "\n";
  opad.strict_sync();
}

