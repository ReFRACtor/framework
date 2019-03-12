#ifndef ABSORBER_VMR_LEVEL_LOG_H
#define ABSORBER_VMR_LEVEL_LOG_H

#include "absorber_vmr_imp_base.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation represents VMR values in log space
  in the state vector.
*******************************************************************/
class AbsorberVmrLevelLog : public AbsorberVmrImpBase {
public:
    AbsorberVmrLevelLog(const boost::shared_ptr<Pressure>& Press,
                        const blitz::Array<double, 1>& Vmr,
                        const blitz::Array<bool, 1>& Vmr_flag,
                        const std::string& Gas_name);

    virtual ~AbsorberVmrLevelLog() = default;

    virtual void print(std::ostream& Os) const;
    
    virtual std::string sub_state_identifier() const
    {
        return "absorber_log_levels/" + gas_name();
    }

    virtual std::string state_vector_name_i(int i) const
    {
        return gas_name() + " Log VMR for Press Lvl " +
               boost::lexical_cast<std::string>(i + 1);
    }

    virtual boost::shared_ptr<AbsorberVmr> clone() const
    {
        return clone(boost::shared_ptr<Pressure>());
    }

    virtual boost::shared_ptr<AbsorberVmr> clone(const boost::shared_ptr<Pressure>& Press) const;

    //-----------------------------------------------------------------------
    /// VMR on the pressure grid. This is just exp(coeff.value), but this is
    /// useful for introspection.
    //-----------------------------------------------------------------------
    blitz::Array<double, 1> vmr_profile() const { return blitz::Array<double, 1>(exp(coeff.value())); }


protected:
    virtual void calc_vmr() const;
};
}
#endif
