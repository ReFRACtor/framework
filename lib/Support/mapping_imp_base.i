// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "mapping_imp_base.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingImpBase);

namespace FullPhysics {

class MappingImpBase : public virtual GenericObject {
public:
  virtual ~MappingImpBase() {};
  virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const;
  virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const = 0;
  virtual std::string name() = 0;
  virtual boost::shared_ptr<MappingImpBase> clone() const = 0;
};
}
