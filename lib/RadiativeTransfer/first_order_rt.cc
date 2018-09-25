#include "first_order_rt.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

FirstOrderRt::FirstOrderRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                           const blitz::Array<double, 1>& Sza, 
                           const blitz::Array<double, 1>& Zen, 
                           const blitz::Array<double, 1>& Azm,
                           bool do_solar, bool do_thermal)
{
    // Watch atmosphere for changes, so we clear cache if needed.
    atm->add_observer(*this);


}

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

blitz::Array<double, 1> FirstOrderRt::stokes_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
}
//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

ArrayAd<double, 1> FirstOrderRt::stokes_and_jacobian_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
}

void FirstOrderRt::print(std::ostream& Os, bool Short_form) const
{
}
