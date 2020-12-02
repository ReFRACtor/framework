%include "fp_common.i"
%{
#include "oss_retrieval_flags.h"
%}

%base_import(generic_object)

%fp_shared_ptr(FullPhysics::OssRetrievalFlags);

namespace FullPhysics {
class OssRetrievalFlags : public GenericObject {
public:
    OssRetrievalFlags(const blitz::Array<int, 1>& Temp_levels, 
                      const blitz::Array<bool, 1>& Skin_temp_flag,
                      const std::vector<blitz::Array<int, 1> >& Gas_levels,
                      const blitz::Array<int, 1>& Emissivity_flags,
                      const blitz::Array<int, 1>& Reflectivity_flags); 

    virtual ~OssRetrievalFlags();

    %python_attribute(temp_levels, blitz::Array<int, 1>);
    %python_attribute(skin_temp_sensors, blitz::Array<bool, 1>);
    %python_attribute(gas_levels, std::vector<blitz::Array<int, 1>>);
    %python_attribute(emissivity_flags, blitz::Array<int, 1>);
    %python_attribute(reflectivity_flags, blitz::Array<int, 1>); 

    %python_attribute(num_total_flags, int);
};
}
