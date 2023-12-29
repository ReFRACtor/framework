// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"

%{
#include "sample_grid.h"
#include "sub_state_vector_array.h"
%}

%fp_shared_ptr(FullPhysics::SampleGrid)
namespace FullPhysics {
    class SampleGrid;
}

%base_import(state_vector_observer)

%import "spectral_domain.i"
%import "sub_state_vector_array.i"

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::SampleGrid>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::SampleGrid>)

namespace FullPhysics {
    
%template(ObservableSampleGrid) FullPhysics::Observable<SampleGrid>;
%template(ObserverSampleGrid) FullPhysics::Observer<SampleGrid>;

%feature("director") SampleGrid;

class SampleGrid: virtual public StateVectorObserver,
public Observable<SampleGrid> {
public:
  virtual ~SampleGrid();
  virtual void add_observer(Observer<SampleGrid>& Obs);
  virtual void remove_observer(Observer<SampleGrid>& Obs);
  virtual boost::shared_ptr<SampleGrid> clone() const = 0;
  // In Python you need to override the renamed function _v_sample_grid
  %python_attribute_abstract(sample_grid, SpectralDomain);
  %python_attribute(pixel_grid, SpectralDomain);
  virtual void state_vector_name(const StateVector& Sv, 
                                 blitz::Array<std::string, 1>& Sv_name) const;
  virtual void notify_update(const StateVector& Observed_object);
  virtual void notify_add(StateVector& Observed_object);
  virtual void notify_remove(StateVector& Observed_object);
  %python_attribute_with_set(sv_name, std::vector<std::string>);
  virtual std::string desc() const;
  std::string print_to_string() const;
  %pickle_serialization();
};
}

%template(vector_sample_grid) std::vector<boost::shared_ptr<FullPhysics::SampleGrid> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(sample_grid, SampleGrid)

// List of things "import *" will include
%python_export("SampleGrid", "ObservableSampleGrid", "ObserverSampleGrid", "vector_sample_grid");
