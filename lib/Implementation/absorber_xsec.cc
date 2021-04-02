#include "absorber_xsec.h"

#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

AbsorberXSec::AbsorberXSec(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
                           const boost::shared_ptr<Pressure>& Press,
                           const boost::shared_ptr<Temperature>& Temp,
                           const std::vector<boost::shared_ptr<Altitude> >& Alt,
                           const std::vector<std::string>& XSec_filenames,
                           const boost::shared_ptr<Constant>& C)
: press(Press), temp(Temp), alt(Alt), vmr(Vmr), xsec_filenames(XSec_filenames), c(C)
{
}

void AbsorberXSec::fill_cache(const AbsorberXSec& T)
{
}

ArrayAd<double, 2> AbsorberXSec::optical_depth_each_layer(double wn, int spec_index) const
{
}

void AbsorberXSec::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "AbsorberXSec:\n";
    for(int i = 0; i < (int) vmr.size(); ++i) {
        Os << "  Absorber VMR[" << i << "]:\n";
        opad << *vmr[i] << "\n";
        opad.strict_sync();
    }
}

boost::shared_ptr<Absorber> AbsorberXSec::clone() const
{
}

boost::shared_ptr<AbsorberVmr> AbsorberXSec::absorber_vmr(const std::string& Gas_name) const
{
    int i = gas_index(Gas_name);
    if(i < 0) {
        Exception err;
        err << "Gas named " << Gas_name << " not present in vmr list";
        throw err;
    }
    return vmr[i];
}
