#include <boost/make_shared.hpp>
#include "absorber_vmr_level_log.h"
#include "mapping_log.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrLevelLog::AbsorberVmrLevelLog(const boost::shared_ptr<Pressure>& Press,
                                         const blitz::Array<double, 1>& Vmr,
                                         const blitz::Array<bool, 1>& Vmr_flag,
                                         const std::string& Gas_name)
: AbsorberVmrLevel(Press, Vmr, Vmr_flag, Gas_name, boost::make_shared<MappingLog>())
{
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevelLog::clone() const
{
    return boost::shared_ptr<AbsorberVmr>(new AbsorberVmrLevelLog(press->clone(), Array<double, 1>(exp(coeff.value())), used_flag, gas_name()));
}