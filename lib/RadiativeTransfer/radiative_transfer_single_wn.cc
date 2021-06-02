#include "radiative_transfer_single_wn.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RadiativeTransferSingleWn::serialize(Archive & ar,
					  const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransferFixedStokesCoefficient)
    & FP_NVP(atm);
}

FP_IMPLEMENT(RadiativeTransferSingleWn);
#endif

Array<double, 2> 
RadiativeTransferSingleWn::stokes(const SpectralDomain& Spec_domain,
				  int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));

  Array<double, 1> wn(Spec_domain.wavenumber());
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);
  Array<double, 2> res(wn.rows(), number_stokes());

  for(int i = 0; i < wn.rows(); ++i) {
    res(i, Range::all()) = stokes_single_wn(wn(i), Spec_index);
    if(disp)	
      *disp += 1;
  }

  Logger::info() << atm->timer_info();

  return res;
}

ArrayAd<double, 2> 
RadiativeTransferSingleWn::stokes_and_jacobian(const SpectralDomain& Spec_domain,
					       int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));

  Array<double, 1> wn(Spec_domain.wavenumber());
  if(wn.rows() < 1)		// Handle degenerate case.
    return ArrayAd<double, 2>(0,number_stokes(),0);
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);
  ArrayAd<double, 1> t =
    stokes_and_jacobian_single_wn(wn(0), Spec_index);
  ArrayAd<double, 2> res(wn.rows(), number_stokes(), t.number_variable());
  res(0, Range::all()) = t;

  // Advance progress bar for first computation used to size the result
  if(disp)
      *disp += 1;

  for(int i = 1; i < wn.rows(); ++i) {
    res(i, Range::all()) = stokes_and_jacobian_single_wn(wn(i), Spec_index);
    if(disp)	
      *disp += 1;
  }

  Logger::info() << atm->timer_info();

  return res;
}

void RadiativeTransferSingleWn::print(std::ostream& Os, bool Short_form) const
{
  RadiativeTransferFixedStokesCoefficient::print(Os, Short_form);
  OstreamPad opad(Os, "  ");
  if(!Short_form) {
    Os << "\nAtmosphere:\n";
    opad << *atm << "\n";
    opad.strict_sync();
  }
}
