// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "spectrum_effect.h"
#include "sub_state_vector_array.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}
%import "observer.i"

%base_import(state_vector_observer)

%import "sub_state_vector_array.i"
%import "spectrum.i"
%import "forward_model_spectral_grid.i"

%nodefaultctor FullPhysics::SubStateVectorArray<SpectrumEffect>;
%fp_shared_ptr(FullPhysics::SpectrumEffect);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::SpectrumEffect>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::SpectrumEffect>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::SpectrumEffect>)

namespace FullPhysics {
  class SpectrumEffect;
}

namespace FullPhysics {
%template(ObservableSpectrumEffect) FullPhysics::Observable<SpectrumEffect>;
%template(ObserverSpectrumEffect) FullPhysics::Observer<SpectrumEffect>;
class SpectrumEffect : virtual public StateVectorObserver,
                               public Observable<SpectrumEffect> {
public:
  virtual ~SpectrumEffect();
  std::string print_to_string() const;
  std::string print_parent() const;
  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const = 0;
  virtual boost::shared_ptr<SpectrumEffect> clone() const = 0;
  %python_attribute_abstract(name, std::string);
  virtual void add_observer(Observer<SpectrumEffect>& Obs);
  virtual void remove_observer(Observer<SpectrumEffect>& Obs);
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization();
};

%template(SubStateVectorArraySpectrumEffect) FullPhysics::SubStateVectorArray<SpectrumEffect>;
}

//%fp_shared_ptr(std::vector<boost::shared_ptr<FullPhysics::SpectrumEffect> >);
//%fp_shared_ptr(std::vector<std::vector<boost::shared_ptr<FullPhysics::SpectrumEffect> > >);

%template(vector_spectrum_effect) std::vector<boost::shared_ptr<FullPhysics::SpectrumEffect> >;
%template(vector_vector_spectrum_effect) std::vector<std::vector<boost::shared_ptr<FullPhysics::SpectrumEffect> > >;
