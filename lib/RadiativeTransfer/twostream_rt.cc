#include "twostream_rt.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void TwostreamRt::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRt);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void TwostreamRt::save(Archive &ar,
		       const unsigned int UNUSED(version)) const
{
  // Save a little extra information if we aren't already saving full
  // state.
  if(!SpurrRt::serialize_full_state) {
    bool do_full_quadrature_ = rt_driver()->do_full_quadrature();
    int number_layer_ = atm->number_layer();
    ar & FP_NVP_(do_full_quadrature) & FP_NVP_(number_layer);
  }
}

template<class Archive>
void TwostreamRt::load(Archive &ar,
		       const unsigned int UNUSED(version)) 
{
  if(!rt_driver_) {
    bool do_full_quadrature_;
    int number_layer_;
    ar & FP_NVP_(do_full_quadrature) & FP_NVP_(number_layer);
    rt_driver_.reset(new TwostreamRtDriver
     (number_layer_, surface_type(), do_full_quadrature_,
      do_solar_sources, do_thermal_emission));
  }
}

FP_IMPLEMENT(TwostreamRt);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TwostreamRt, RadiativeTransfer)
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&>())
.def(luabind::constructor<const boost::shared_ptr<RtAtmosphere>&,
                          const boost::shared_ptr<StokesCoefficient>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&, 
                          const blitz::Array<double, 1>&,
                          bool>())
REGISTER_LUA_END()
#endif

// For debugging purposes, it can be useful to dump out the input
// used by this class. We'll leave this code in place, in case we
// need it again, but normally this turned off.
const bool dump_data = false;

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm Atmosphere class with optical properties to model
/// \param Stokes_coef Multiplier for each stokes component
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param do_fullquadrature false only for comparison against LIDORT 
//-----------------------------------------------------------------------

TwostreamRt::TwostreamRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                         const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                         const blitz::Array<double, 1>& Sza,
                         const blitz::Array<double, 1>& Zen,
                         const blitz::Array<double, 1>& Azm, 
                         bool do_fullquadrature, bool do_solar, bool do_thermal)
: SpurrRt(Atm, Stokes_coef, Sza, Zen, Azm, do_solar, do_thermal)
{   
  rt_driver_.reset(new TwostreamRtDriver(atm->number_layer(), surface_type(), do_fullquadrature, do_solar_sources, do_thermal_emission));
  if(dump_data)
    std::cout << "# Nlayer:\n" << atm->number_layer() << "\n"
              << "# Surface type:\n" << surface_type() << "\n"
              << "# do_fullquadrature:\n" << do_fullquadrature << "\n";
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void TwostreamRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "TwostreamRt\n";
  OstreamPad opad1(Os, "  ");
  SpurrRt::print(opad1, Short_form);
  opad1.strict_sync();
  Os << "do_full_quadrature = " << (rt_driver()->do_full_quadrature() ? "true" : "false") << "\n";
  opad1.strict_sync();
}
