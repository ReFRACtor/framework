#ifndef CHAPMAN_BOA_RT_H
#define CHAPMAN_BOA_RT_H
#include "radiative_transfer.h"
#include "atmosphere_standard.h"
#include "chapman_boa.h"
#include "spectral_bound.h"
#include "calculation_cache.h"

namespace FullPhysics {
class ChapmanBoaRT;

class ChapmanBoaRTCache : public CalculationCacheIndexed<ChapmanBoaRT>
{
public:
  ChapmanBoaRTCache() {}
  virtual ~ChapmanBoaRTCache() {}
  void fill_cache(const ChapmanBoaRT& C, int Spec_index);
  boost::shared_ptr<ChapmanBOA> chapman_boa;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
  
/****************************************************************//**
*******************************************************************/

class ChapmanBoaRT : public RadiativeTransfer {
public:
  ChapmanBoaRT(const boost::shared_ptr<AtmosphereStandard>& Atm,
  	       const blitz::Array<double, 1>& Sza);

  ChapmanBoaRT(const boost::shared_ptr<AtmosphereStandard>& Atm,
	       const blitz::Array<double, 1>& Sza, 
	       const SpectralBound& Spec_bound);

  virtual ~ChapmanBoaRT() {}

  virtual int number_stokes() const { return 1; }

  virtual int number_spectrometer() const
  { return sza.extent(blitz::firstDim);}

  //-----------------------------------------------------------------------
  /// Pointer to the Atmosphere class we are using.
  //-----------------------------------------------------------------------

  const boost::shared_ptr<AtmosphereStandard>& atmosphere_ptr() const {return atm;}

  // See description in base class
  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const;
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian(const SpectralDomain& Spec_domain, int Spec_index) const;
 
  virtual void print(std::ostream& Os, bool Short_form = false) const;
  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    res.push_back(atm);
    return res;
  }
private:
  friend ChapmanBoaRTCache;
  mutable std::vector<ChapmanBoaRTCache> cache;

  /// Window range for use when picking refraction method
  SpectralBound spec_bound;

  /// Atmosphere object we are using.
  boost::shared_ptr<AtmosphereStandard> atm;

  /// Solar zenith angles per spectrometer
  blitz::Array<double, 1> sza;
  ChapmanBoaRT() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ChapmanBoaRT);
#endif
