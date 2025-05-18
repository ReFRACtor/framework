%include "fp_common.i"

%{
#include "state_mapping_interpolate.h"
#include "linear_interpolate.h"
#include "log_interpolate.h"
%}

%base_import(state_mapping)

%import "array_ad.i"
%import "pressure.i"

namespace FullPhysics {
    template <class Interp> class StateMappingInterpolate;
}

%define %map_interpolate_template(NAME, INTERP)

%fp_shared_ptr(FullPhysics::StateMappingInterpolate<INTERP>);

namespace FullPhysics {

template <>
class StateMappingInterpolate<INTERP> : public StateMapping {
public:
    virtual ~StateMappingInterpolate() {}
    StateMappingInterpolate(const boost::shared_ptr<Pressure>& PressTo,
                            const boost::shared_ptr<Pressure>& PressFrom);
    virtual boost::shared_ptr<StateMapping> clone() const;
    %pickle_serialization();
};

}

%template(NAME) FullPhysics::StateMappingInterpolate<INTERP>;

%enddef

%map_interpolate_template(StateMappingInterpolateLinearLinear, FullPhysics::LinearInterpolate);
%map_interpolate_template(StateMappingInterpolateLogLinear, FullPhysics::LogLinearInterpolate);
%map_interpolate_template(StateMappingInterpolateLogLog, FullPhysics::LogLogInterpolate);
%map_interpolate_template(StateMappingInterpolateLinearLog, FullPhysics::LinearLogInterpolate);
