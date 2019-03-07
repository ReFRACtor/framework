#include <boost/bind.hpp>
#include "absorber_vmr_level_log.h"
#include "ostream_pad.h"
#include "linear_interpolate.h"
using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrLevelLog::AbsorberVmrLevelLog(const boost::shared_ptr<Pressure>& Press,
                                         const blitz::Array<double, 1>& Vmr,
                                         const blitz::Array<bool, 1>& Vmr_flag,
                                         const std::string& Gas_name)
: AbsorberVmrImpBase(Gas_name, Array<double, 1> (log(Vmr)), Vmr_flag, Press, false)
{
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevelLog::clone(const boost::shared_ptr<Pressure>& Press) const
{
    return boost::shared_ptr<AbsorberVmr>(new AbsorberVmrLevelLog(Press, Array<double, 1>(exp(coeff.value())), used_flag, gas_name()));
}

void AbsorberVmrLevelLog::calc_vmr() const
{
    std::vector<AutoDerivative<double> > plist;
    std::vector<AutoDerivative<double> > vmrlist;

    for(int i = 0; i < press->pressure_grid().rows(); ++i) {
        vmrlist.push_back(exp(coeff(i)));
        plist.push_back(press->pressure_grid()(i).value);
    }

    typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> > lin_type;

    boost::shared_ptr<lin_type> lin (new lin_type(plist.begin(), plist.end(), vmrlist.begin()));

    vmr = boost::bind(&lin_type::operator(), lin, _1);
}

void AbsorberVmrLevelLog::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "AbsorberVmrLevelLog:\n"
       << "  Gas name:       " << gas_name() << "\n"
       << "  Coefficient:\n";
    opad << coeff.value() << "\n";
    opad.strict_sync();
    Os << "  Retrieval Flag:\n";
    opad << used_flag << "\n";
    opad.strict_sync();
}
