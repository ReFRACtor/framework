#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "printable.h"
#include "spectral_window.h"
#include "array_ad.h"
#include "state_vector_observer.h"
#include "spectrum.h"
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This applies a instrument model to radiances.
*******************************************************************/

class Instrument : virtual public StateVectorObserver, 
		   public Observable<Instrument> {
public:
  virtual ~Instrument() {}

  virtual void add_observer(Observer<Instrument>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Instrument>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Clone an Instrument object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Instrument>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Instrument object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Instrument> clone() const = 0;

//-----------------------------------------------------------------------
/// Give number of spectrometers.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// Apply the instrument model to both the radiance and derivatives.
///
/// \param High_resolution_spectrum High resolution spectrum.
/// \param Pixel_list List of pixels to include in radiance
/// \param Spec_index Spectral index
/// \return Spectrum with instrument model applied.
//-----------------------------------------------------------------------

  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// This is the pixel wavenumber/wavelength for each pixel.
//-----------------------------------------------------------------------

  virtual SpectralDomain pixel_spectral_domain(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// This is the amount of spectral points added on each edge of the
/// monochromatic grid to allow for enough points for ILS convolution
//-----------------------------------------------------------------------

  virtual DoubleWithUnit high_res_extension(int Spec_index) const = 0;
  virtual void high_res_extension(int Spec_index, DoubleWithUnit& extension) = 0;
  virtual void print(std::ostream& Os) const {Os << "Instrument";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(Instrument);
FP_EXPORT_OBSERVER_KEY(Instrument);
#endif
