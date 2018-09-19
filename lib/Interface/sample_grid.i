// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"

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
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::SampleGrid>)
%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::SampleGrid>;

namespace FullPhysics {
    
    %template(ObservableSampleGrid) FullPhysics::Observable<SampleGrid>;
    %template(ObserverSampleGrid) FullPhysics::Observer<SampleGrid>;
    
    class SampleGrid: virtual public StateVectorObserver,
    public Observable<SampleGrid> {
    public:
        virtual ~SampleGrid();
        virtual void add_observer(Observer<SampleGrid>& Obs);
        virtual void remove_observer(Observer<SampleGrid>& Obs);
        virtual boost::shared_ptr<SampleGrid> clone() const = 0;
        %python_attribute(sample_grid, virtual SpectralDomain);
        %python_attribute(pixel_grid, virtual SpectralDomain);
    };
    
    %template(SubStateVectorArraySampleGrid) FullPhysics::SubStateVectorArray<FullPhysics::SampleGrid>;
}

%template(vector_sample_grid) std::vector<boost::shared_ptr<FullPhysics::SampleGrid> >;
