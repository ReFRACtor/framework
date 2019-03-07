%include "common.i"

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

    virtual ~AbsorberVmrLevelLog();
    virtual void print(std::ostream& Os) const;
    virtual std::string sub_state_identifier() const;
    virtual std::string state_vector_name_i(int i) const;
    virtual boost::shared_ptr<AbsorberVmr> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> clone(const boost::shared_ptr<Pressure>& Press) const;
};
}
