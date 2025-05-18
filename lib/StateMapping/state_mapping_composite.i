%include "fp_common.i"

%{
#include "state_mapping_composite.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::StateMappingComposite);

namespace FullPhysics {

class StateMappingComposite : public StateMapping  {
public:
    virtual ~StateMappingComposite() {};
    StateMappingComposite(const std::vector<boost::shared_ptr<StateMapping> >& Mappings);
    virtual boost::shared_ptr<StateMapping> clone() const;
    %pickle_serialization();
};
}
