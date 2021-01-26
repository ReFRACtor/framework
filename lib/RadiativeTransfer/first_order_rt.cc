#include "first_order_rt.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void FirstOrderRt::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRt);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void FirstOrderRt::save(Archive &ar,
		   const unsigned int UNUSED(version)) const
{
  // Save a little extra information if we aren't already saving full
  // state.
  if(!SpurrRt::serialize_full_state) {
    int number_stream_ = number_stream();
    int number_moment_ = number_moment();
    ar & FP_NVP_(number_stream) & FP_NVP_(number_moment);
  }
}

template<class Archive>
void FirstOrderRt::load(Archive &ar,
			const unsigned int UNUSED(version)) 
{
  if(!rt_driver_) {
    int number_stream_;
    int number_moment_;
    ar & FP_NVP_(number_stream) & FP_NVP_(number_moment);
    rt_driver_.reset(new FirstOrderDriver
     (atm->number_layer(), surface_type(), number_stream_, number_moment_,
      do_solar_sources, do_thermal_emission));
  }
}

FP_IMPLEMENT(FirstOrderRt);
#endif

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Atm The atmosphere to use
/// \param Stokes_coef The Stokes coefficients.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Number_streams Number of streams to use
/// \param Number_moments Number of moments to use in phase function (-1 to use
///     full resolution
//-----------------------------------------------------------------------

FirstOrderRt::FirstOrderRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                           const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                           const blitz::Array<double, 1>& Sza, 
                           const blitz::Array<double, 1>& Zen, 
                           const blitz::Array<double, 1>& Azm,
                           int Number_streams,
                           int Number_moments,
                           bool do_solar,
                           bool do_thermal)
: SpurrRt(Atm, Stokes_coef, Sza, Zen, Azm, do_solar, do_thermal)
{   
  rt_driver_.reset(new FirstOrderDriver(atm->number_layer(), surface_type(), Number_streams, Number_moments, do_solar_sources, do_thermal_emission));
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void FirstOrderRt::print(std::ostream& Os, bool Short_form) const
{
  Os << "FirstOrderRt\n";
  OstreamPad opad1(Os, "  ");
  SpurrRt::print(opad1, Short_form);
  opad1.strict_sync();
  OstreamPad opad(Os, "    ");
  Os << "  Number stream: " << number_stream() << "\n"
     << "  Number moments: " << number_moment() << "\n";
  opad << "\n";
  opad.strict_sync();
}

