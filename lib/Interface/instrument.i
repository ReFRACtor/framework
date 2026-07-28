%include "fp_common.i"

%{
#include "instrument.h"
%}

%import "observer.i"
%import "spectrum.i"
%import "spectral_domain.i"
%import "double_with_unit.i"
%base_import(state_vector_observer)

%fp_shared_ptr(FullPhysics::Instrument);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Instrument>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Instrument>)

namespace FullPhysics {
  class Instrument;
}

namespace FullPhysics {
%template(ObservableInstrument) FullPhysics::Observable<Instrument>;
%template(ObserverInstrument) FullPhysics::Observer<Instrument>;

// Allow this class to be derived from in Python.
%feature("director") Instrument;

class Instrument : virtual public StateVectorObserver,
		   public Observable<Instrument> {
public:
  virtual ~Instrument();
  std::string print_to_string() const;
  std::string print_parent() const;
  virtual void add_observer(Observer<Instrument>& Obs);
  virtual void remove_observer(Observer<Instrument>& Obs);
  virtual boost::shared_ptr<Instrument> clone() const = 0;
  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const = 0;

  // When implementing directors can not use %python_attribute here or else there will
  // be missing symbol problems in the director
  virtual int number_spectrometer() const = 0;
  virtual SpectralDomain pixel_spectral_domain(int Spec_index) const = 0;
  virtual DoubleWithUnit high_res_extension(int Spec_index) const = 0;
  virtual void high_res_extension(int Spec_index, DoubleWithUnit& extension) = 0;
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(instrument, Instrument)

// List of things "import *" will include
%python_export("Instrument", "ObservableInstrument", "ObserverInstrument");
