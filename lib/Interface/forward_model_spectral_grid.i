// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "forward_model_spectral_grid.h"
#include "instrument.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"
%import "instrument.i"
%import "spectral_window.i"
%import "spectrum_sampling.i"

%fp_shared_ptr(FullPhysics::ForwardModelSpectralGrid);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") ForwardModelSpectralGrid;

class ForwardModelSpectralGrid  : public GenericObject {
public:
  ForwardModelSpectralGrid(
   const boost::shared_ptr<Instrument>& Inst,
   const boost::shared_ptr<SpectralWindow>& Spectral_window,
   const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling);
  ForwardModelSpectralGrid() {}
  std::string print_to_string() const;
  std::string print_parent() const;
  %python_attribute(number_spectrometer, virtual int);
  SpectralDomain low_resolution_grid(int Spec_index) const;
  SpectralDomain high_resolution_grid(int Spec_index) const;
  SpectralDomain high_resolution_interpolated_grid(int Spec_index) const;
  Spectrum interpolate_spectrum(const Spectrum& Spec_in, int Spec_index) const;
  const std::vector<int> pixel_list(int Spec_index) const;
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(forward_model_spectral_grid, ForwardModelSpectralGrid)

// List of things "import *" will include
%python_export("ForwardModelSpectralGrid");
