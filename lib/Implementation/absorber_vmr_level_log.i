%include "fp_common.i"

%{
#include "absorber_vmr_level_log.h"
%}

%base_import(absorber_vmr_imp_base)
%import "pressure.i"

%fp_shared_ptr(FullPhysics::AbsorberVmrLevelLog)

namespace FullPhysics {
class AbsorberVmrLevelLog : public AbsorberVmrImpBase {
public:
    AbsorberVmrLevelLog(const boost::shared_ptr<Pressure>& Press,
                        const blitz::Array<double, 1>& Vmr,
                        const blitz::Array<bool, 1>& Vmr_flag,
                        const std::string& Gas_name);
    virtual boost::shared_ptr<AbsorberVmr> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> clone(const boost::shared_ptr<Pressure>& Press) const;
    %python_attribute(sub_state_identifier, std::string);
    virtual std::string state_vector_name_i(int i) const;
    %python_attribute(vmr_profile, blitz::Array<double, 1>);
protected:
    virtual void calc_vmr() const;
};
}
