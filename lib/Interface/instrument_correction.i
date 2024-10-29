// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

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

%nodefaultctor FullPhysics::SubStateVectorArray<InstrumentCorrection>;
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
class InstrumentCorrection : virtual public StateVectorObserver,
			     public Observable<InstrumentCorrection> {
public:
  virtual ~InstrumentCorrection();
  std::string print_to_string() const;
  std::string print_parent() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const = 0;
  virtual boost::shared_ptr<InstrumentCorrection> clone() const = 0;
  %pickle_serialization();
};

%template(SubStateVectorArrayInstrumentCorrection) FullPhysics::SubStateVectorArray<InstrumentCorrection>;
}

%template(vector_instrument_correction) std::vector<boost::shared_ptr<FullPhysics::InstrumentCorrection> >;
%template(vector_vector_instrument_correction) std::vector<std::vector<boost::shared_ptr<FullPhysics::InstrumentCorrection> > >;
