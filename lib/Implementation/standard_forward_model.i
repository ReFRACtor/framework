// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "standard_forward_model.h"
#include "ils_instrument.h"
#include "pressure.h"
#include "spectrum.h"

// Needed for type conversions in SWIG
#include "sub_state_vector_array.h"
%}
%base_import(forward_model)
%base_import(named_spectrum)
%import "forward_model_spectral_grid.i"
%import "instrument.i"
%import "radiative_transfer.i"
%import "spectral_window.i"
%import "spectrum.i"
%import "spectrum_sampling.i"
%import "spectrum_effect.i"

%fp_shared_ptr(FullPhysics::StandardForwardModel);

namespace FullPhysics {

// Allow these classes to be derived from in Python.
%feature("director") StandardForwardModel;

// Note, a class that is derived from in python needs to declare every virtual function that
// can be called on it, even if all that happens is the base class
// to a director is called. This is because this class is used to
// create the SwigDirector class, and this class needs each of the member functions to
// direct things properly. It is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
  
class StandardForwardModel : public ForwardModel,
   public Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> > {
public:
  StandardForwardModel(
   const boost::shared_ptr<Instrument>& Inst,
   const boost::shared_ptr<SpectralWindow>& Spectral_window,
   const boost::shared_ptr<RadiativeTransfer>& Rt,
   const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling,
   const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& Spectrum_effect = 
		  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >());
  virtual Spectrum radiance(int sensor_index, bool Skip_jacobian = false) 
    const;
  virtual Spectrum radiance_all(bool skip_jacobian = false) const;
  %python_attribute(num_channels, virtual int)
  virtual SpectralDomain spectral_domain(int sensor_index) const;
  virtual void setup_grid();
  virtual SpectralDomain::TypePreference spectral_domain_type_preference() const;
  %python_attribute_with_set(instrument, boost::shared_ptr<Instrument>);
  %python_attribute_with_set(spectral_window, boost::shared_ptr<SpectralWindow>);
  %python_attribute_with_set(radiative_transfer, boost::shared_ptr<RadiativeTransfer>);
  %python_attribute_with_set(spectrum_sampling, boost::shared_ptr<SpectrumSampling>);
  %python_attribute_with_set(spectral_grid, boost::shared_ptr<ForwardModelSpectralGrid>);
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);

  Spectrum apply_spectrum_corrections(const Spectrum& highres_spec, int sensor_index) const;

  virtual void add_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs); 
  virtual void remove_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs);

  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const;

  std::string print_to_string() const;
  virtual std::string desc() const;

  // vector of vector for SpectrumEffect is kind of a pain in
  // python. So just brute force a conversion.
  %extend {
    int _speceff_size() { return (int) $self->spectrum_effect().size(); }
    int speceff_size2(int i) { return (int) $self->spectrum_effect()[i].size();}
    boost::shared_ptr<SpectrumEffect> speceff_val(int i, int j)
    {return $self->spectrum_effect()[i][j];}
  }
  %pythoncode {
@property
def speceff_size(self):
  return self._speceff_size()

@property
def spectrum_effect(self):
  res = []
  for i in range(self.speceff_size):
     res2 = []
     for j in range(self.speceff_size2(i)):
        res2.append(self.speceff_val(i,j))
     res.append(res2)
  return res
  }
  %pickle_serialization();
protected:
  StandardForwardModel();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(standard_forward_model, StandardForwardModel)

// List of things "import *" will include
%python_export("StandardForwardModel");
