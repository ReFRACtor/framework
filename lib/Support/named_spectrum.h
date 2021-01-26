#ifndef NAMED_SPECTRUM_H
#define NAMED_SPECTRUM_H

#include "spectrum.h"
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
 Adds name and spec index fields to a Spectrum. Useful for sending
 Spectrum files to output files.
*******************************************************************/
class NamedSpectrum: public Spectrum {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

  NamedSpectrum(const SpectralDomain& Spec_domain, 
                const SpectralRange& Spec_range, const std::string& Name,
                int Index)
    : Spectrum(Spec_domain, Spec_range), name_(Name), index_(Index) {}

  NamedSpectrum(const Spectrum& Spec, const std::string& Name, int Index)
    : Spectrum(Spec.spectral_domain(), Spec.spectral_range()),
      name_(Name), index_(Index) {}

//-----------------------------------------------------------------------
/// Name that makes this a named spectrum
//-----------------------------------------------------------------------

  virtual const std::string& name() const {return name_;}

//-----------------------------------------------------------------------
/// An reference index for the spectrum, ie a spectrometer index 
//-----------------------------------------------------------------------

  virtual int index() const {return index_;}

  void print(std::ostream& Os) 
  { 
    Os << "NamedSpectrum:" << std::endl
       << "  name: " << name_ << std::endl
       << "  index: " << index_;
  }

  /// Default constructor needed for SWIG
  NamedSpectrum() {}

private:
  std::string name_;
  int index_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<NamedSpectrum> NamedSpectrumPtr;
typedef std::vector<boost::shared_ptr<NamedSpectrum> > NamedSpectrumPtrVec;
  
}

FP_EXPORT_KEY(NamedSpectrum);
FP_EXPORT_OBSERVER_KEY(NamedSpectrumPtr);
FP_EXPORT_OBSERVER_KEY(NamedSpectrumPtrVec);
#endif
