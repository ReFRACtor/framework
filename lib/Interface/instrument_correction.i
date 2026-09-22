%include <std_vector.i>
%include "fp_common.i"

%{
#include "instrument_correction.h"
#include "sub_state_vector_array.h"
%}

%base_import(state_vector_observer)

%import "spectral_domain.i"
%import "spectral_range.i"
%import "sub_state_vector_array.i"

%template(vector_instrument_correction) std::vector<boost::shared_ptr<FullPhysics::InstrumentCorrection> >;
%template(vector_vector_instrument_correction) std::vector<std::vector<boost::shared_ptr<FullPhysics::InstrumentCorrection> > >;
%fp_shared_ptr(FullPhysics::InstrumentCorrection);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::InstrumentCorrection>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::InstrumentCorrection>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::InstrumentCorrection>)

namespace FullPhysics {
  class InstrumentCorrection;
}

namespace FullPhysics {
%template(ObservableInstrumentCorrection) FullPhysics::Observable<InstrumentCorrection>;
%template(ObserverInstrumentCorrection) FullPhysics::Observer<InstrumentCorrection>;

// Allow this class to be derived from in Python.
%feature("director") InstrumentCorrection;

class InstrumentCorrection : virtual public StateVectorObserver,
			     public Observable<InstrumentCorrection> {
public:
  virtual ~InstrumentCorrection();
  virtual void add_observer(Observer<InstrumentCorrection>& Obs);
  virtual void remove_observer(Observer<InstrumentCorrection>& Obs);
  std::string print_to_string() const;
  std::string print_parent() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const = 0;
  virtual boost::shared_ptr<InstrumentCorrection> clone() const = 0;
  %pickle_serialization();
};

// Allow the SubStateVectorArray<InstrumentCorrection> base to be derived
// from in Python as well. As with other director classes derived from
// SubStateVectorArray, every virtual function that can be called on it needs
// to be listed explicitly (see %sub_state_virtual_func).
%feature("director") SubStateVectorArray<InstrumentCorrection>;

%template(SubStateVectorArrayInstrumentCorrection) SubStateVectorArray<InstrumentCorrection>;
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(instrument_correction, InstrumentCorrection)

// List of things "import *" will include
%python_export("InstrumentCorrection", "ObservableInstrumentCorrection", "ObserverInstrumentCorrection", "SubStateVectorArrayInstrumentCorrection");

