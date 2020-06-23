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

class Mapping : public GenericObject {
public:
  virtual ~Mapping();
  virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff) const;
  virtual ArrayAd<double, 1> retrieval_init
  (const ArrayAd<double, 1>& initial_coeff) const;
  %python_attribute(name, std::string);
  virtual boost::shared_ptr<Mapping> clone() = 0;
  std::string print_to_string();
  %pickle_serialization();
};
}
