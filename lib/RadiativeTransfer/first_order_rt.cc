#include "first_order_rt.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void FirstOrderRt::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrRt);
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
  rt_driver_.reset(new FirstOrderDriver(atm->number_layer(), surface_type(), Number_streams, Number_moments, do_solar, do_thermal));
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

