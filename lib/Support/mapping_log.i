%include "fp_common.i"

%{
#include "mapping_log.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingLog);


namespace FullPhysics {

%feature("notabstract") MappingLog;

class MappingLog : public Mapping {
public:
    MappingLog();
    virtual ~MappingLog();
    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const;
    virtual std::string name() const;
    virtual boost::shared_ptr<Mapping> clone() const;
};
}
