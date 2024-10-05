#ifndef SPECTRUM_SAMPLING_H
#define SPECTRUM_SAMPLING_H

#include "printable.h"
#include "spectral_domain.h"
#include "spectrum.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This determines the sampling of the monochromatic spectrum
  that should be used for each of the spectrum indexes.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrumdoxygen for a description of each
  of these.
*******************************************************************/
class SpectrumSampling : public Printable<SpectrumSampling> {
public:
  virtual ~SpectrumSampling() {}

//-----------------------------------------------------------------------
/// Number of spectrometers we have.
//-----------------------------------------------------------------------

  int number_spectrometer() const { return nspectrometer; }

//-----------------------------------------------------------------------
/// Wavenumbers/Wavelengths to use for the given spectrometer. We pass
/// in the low resolution grid that we are going to generate after the
/// ILS convolution, along with the edge extension amount so we can generate
/// the high resolution points needed to supply the ILS.
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain(int spec_index, 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Edge_extension) const = 0;

//-----------------------------------------------------------------------
/// The interpolated spectral domain. The default is that this is just
/// the same as spectral_domain, but derived classes can supply a
/// different implementation if it is doing nonuniform sampling.
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain_interpolated(int Spec_index, 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Edge_extension) const
  { return spectral_domain(Spec_index, Lowres_grid, Edge_extension); }

//-----------------------------------------------------------------------
/// Indicate if spectral_domain and spectral_domain_interpolated are
/// different at all.
//-----------------------------------------------------------------------

  virtual bool need_interpolation(int UNUSED(Spec_index)) const { return false; }

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "SpectrumSampling";}

//-----------------------------------------------------------------------
/// We have some fairly nested object hierarchies. It can be useful to
/// be able to search this for things (e.g., which Pressure object is
/// used by a ForwardModel?). This returns a list of subobjects
/// "owned" by this object.
//-----------------------------------------------------------------------

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    return res;
  }
protected:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
  SpectrumSampling(int num_spectrometer) : nspectrometer(num_spectrometer) {}

//-----------------------------------------------------------------------
/// Default constructor, derived classes should set nspectrometer.
//-----------------------------------------------------------------------
  SpectrumSampling() {}

  int nspectrometer;		//< Derived classes should set this.
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Handle the degenerate case where we want the high and low resolution
  to be the same.
*******************************************************************/
  
class IdentitySpectrumSampling: public SpectrumSampling {
public:
  IdentitySpectrumSampling(int nspec) : SpectrumSampling(nspec) {}
  virtual ~IdentitySpectrumSampling() {}
  virtual SpectralDomain spectral_domain(int UNUSED(spec_index), 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& UNUSED(Edge_extension)) const
  {
    return Lowres_grid;
  }
  virtual void print(std::ostream& Os) const {Os << "IdentitySpectrumSampling";}
private:
  IdentitySpectrumSampling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
    
}

FP_EXPORT_KEY(SpectrumSampling);
FP_EXPORT_KEY(IdentitySpectrumSampling);
#endif
