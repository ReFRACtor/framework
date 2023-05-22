// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "measured_radiance_field.h"
%}

%base_import(observer)
%base_import(state_vector)
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::MeasuredRadianceField)
namespace FullPhysics {
  class MeasuredRadianceField;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::MeasuredRadianceField>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::MeasuredRadianceField>)

namespace FullPhysics {
%template(ObservableMeasuredRadianceField) FullPhysics::Observable<FullPhysics::MeasuredRadianceField>;
%template(ObserverMeasuredRadianceField) FullPhysics::Observer<FullPhysics::MeasuredRadianceField>;
  
class MeasuredRadianceField : public Observable<MeasuredRadianceField> {
public:
  virtual ~MeasuredRadianceField();
  virtual void add_observer(Observer<MeasuredRadianceField>& Obs); 
  virtual void remove_observer(Observer<MeasuredRadianceField>& Obs);
  std::string print_to_string() const;
  virtual boost::shared_ptr<Spectrum> measured_radiance_field(int spec_index) const = 0;
  virtual boost::shared_ptr<MeasuredRadianceField> clone() const;
  virtual void print(std::ostream& Os) const;
  %pickle_serialization();
};
}

