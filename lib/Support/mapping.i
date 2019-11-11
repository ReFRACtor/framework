// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "mapping.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::Mapping);

namespace FullPhysics {

class Mapping : public virtual GenericObject {
public:
  virtual ~Mapping() {};
  virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const = 0;
  virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const = 0;
  virtual std::string name() = 0;
  virtual boost::shared_ptr<Mapping> clone() const = 0;
};
}
