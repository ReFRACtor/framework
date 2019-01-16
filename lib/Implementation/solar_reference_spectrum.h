#ifndef SOLAR_REF_SPEC_H
#define SOLAR_REF_SPEC_H

#include <boost/shared_ptr.hpp>

#include "solar_model.h"
#include "solar_doppler_shift.h"

#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This applies a solar model to radiances to model the incoming solar
  irradiance.

  This implementation uses reference solar spectrum data that 
  includes continuum and absorption effects together.

  Additionally a SolarDopplerShift object applies a correction for
  the distance from the sun.
*******************************************************************/

class SolarReferenceSpectrum : public SolarModel {
public:

    //-----------------------------------------------------------------------
    /// Uses the supplied Spectrum object as the reference spectrum along
    /// along with an optional SolarDopplerShift object
    //-----------------------------------------------------------------------

    SolarReferenceSpectrum(const boost::shared_ptr<Spectrum>& reference_spectrum,
                           const boost::shared_ptr<SolarDopplerShift>& doppler_shift = NULL);

    virtual ~SolarReferenceSpectrum() = default;

    //-----------------------------------------------------------------------
    /// Clone a SolarReferenceSpectrum object
    //----------------------------------------------------------------------- 

    virtual boost::shared_ptr<SpectrumEffect> clone() const {
        boost::shared_ptr<SpectrumEffect> res (new SolarReferenceSpectrum(ref_spec_orig, doppler_shift_));
        return res;
    }
 
    //-----------------------------------------------------------------------
    /// The SolarDopplerShift object used by this class.
    //-----------------------------------------------------------------------

    virtual const boost::shared_ptr<SolarDopplerShift>& doppler_shift() const { return doppler_shift_; }

    //-----------------------------------------------------------------------
    /// The SolarContinuumSpectrum object used by this class, as a ptr.
    //-----------------------------------------------------------------------

    virtual void print(std::ostream& Os) const;
    virtual Spectrum solar_spectrum(const SpectralDomain& spec_domain) const;

private:

    boost::shared_ptr<Spectrum> ref_spec_orig; // Stored for cloning
    LinearInterpolate<double, double> ref_spec_interp;
    boost::shared_ptr<SolarDopplerShift> doppler_shift_;

};
}
#endif
