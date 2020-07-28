#ifndef STANDARD_FORWARD_MODEL_H
#define STANDARD_FORWARD_MODEL_H

#include "forward_model.h"
#include "forward_model_spectral_grid.h"
#include "radiative_transfer.h"
#include "spectrum_effect.h"
#include "named_spectrum.h"
#include "spectral_window.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "observer.h"
#include <vector>

namespace FullPhysics {

/****************************************************************//**
  This is the forward model used form GOSAT/OCO. This is fairly
  general, we may want to rename this at some point if we use it for
  other instruments.
*******************************************************************/

class StandardForwardModel : public ForwardModel,
			     public Observable<boost::shared_ptr<NamedSpectrum> > {
public:
  StandardForwardModel(
		       const boost::shared_ptr<Instrument>& Inst,
		       const boost::shared_ptr<SpectralWindow>& Spectral_window,
		       const boost::shared_ptr<RadiativeTransfer>& Rt,
		       const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling,
		       const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& Spectrum_effect =
		       std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >());

  virtual ~StandardForwardModel() {}

  virtual void setup_grid()
  {
    g.reset(new ForwardModelSpectralGrid(inst, swin, spectrum_sampling_));
  }

  virtual int num_channels() const
  {
    return swin->number_spectrometer();
  }

  virtual SpectralDomain spectral_domain(int channel_index) const
  {
    if(!g) {
      throw Exception ("setup_grid needs to be called before calling spectral_domain");
    }

    return g->low_resolution_grid(channel_index);
  }

  virtual SpectralDomain::TypePreference spectral_domain_type_preference() const
  {
    return inst->pixel_spectral_domain(0).type_preference();
  }

  virtual Spectrum radiance(int channel_index, bool Skip_jacobian = false) const;

  virtual void print(std::ostream& Os) const;

  // These are useful for python, so make available here.
  const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& spectrum_effect() const
  {
    return spec_effect;
  }

  const boost::shared_ptr<Instrument>& instrument() const
  {
    return inst;
  }

  void instrument(const boost::shared_ptr<Instrument>& V)
  {
    inst = V;
  }

  const boost::shared_ptr<SpectralWindow>& spectral_window() const
  {
    return swin;
  }

  void spectral_window(const boost::shared_ptr<SpectralWindow>& V)
  {
    swin = V;
  }

  const boost::shared_ptr<RadiativeTransfer>& radiative_transfer() const
  {
    return rt;
  }
  void radiative_transfer(const boost::shared_ptr<RadiativeTransfer>& V)
  {
    rt = V;
  }

  const boost::shared_ptr<SpectrumSampling>& spectrum_sampling() const
  {
    return spectrum_sampling_;
  }

  void spectrum_sampling(const boost::shared_ptr<SpectrumSampling>& V)
  {
    spectrum_sampling_ = V;
  }

  Spectrum apply_spectrum_corrections(const Spectrum& highres_spec, int channel_index) const;

  const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grid() const
  {
    return g;
  }

  /// Required observable functions
  virtual void add_observer(Observer<boost::shared_ptr<NamedSpectrum> > & Obs)
  {
    add_observer_do(Obs);
  }

  virtual void remove_observer(Observer<boost::shared_ptr<NamedSpectrum> >& Obs)
  {
    remove_observer_do(Obs);
  }

  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int channel_index) const;

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const;
private:
  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > spec_effect;
  boost::shared_ptr<Instrument> inst;
  boost::shared_ptr<SpectralWindow> swin;
  boost::shared_ptr<RadiativeTransfer> rt;
  boost::shared_ptr<SpectrumSampling> spectrum_sampling_;
  boost::shared_ptr<ForwardModelSpectralGrid> g;
  StandardForwardModel() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StandardForwardModel);
#endif
