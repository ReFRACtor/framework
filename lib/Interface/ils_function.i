%include "fp_common.i"

%{
#include "ils_function.h"
%}
%base_import(generic_object)

%import "auto_derivative.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::IlsFunction);

namespace FullPhysics {

%feature("director") IlsFunction;

class IlsFunction : public GenericObject {
public:
  virtual ~IlsFunction();
  virtual std::string desc() const;
  std::string print_to_string() const;
  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& OUTPUT) const = 0;

  // Needed for directors, can not use %python_attribute here or else there will
  // be missing symbol problems in the director
  virtual std::string band_name() const = 0;
  virtual std::string hdf_band_name() const;
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%{
// Needed by code below, can't easily figure these names out
// automatically so just include here
#include "ils_function_wrap.h"
%}
%fp_director_serialization(IlsFunction)

// List of things "import *" will include
%python_export("IlsFunction");
