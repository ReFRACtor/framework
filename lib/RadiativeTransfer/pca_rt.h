#ifndef PCA_RT_H
#define PCA_RT_H

#include "radiative_transfer_fixed_stokes_coefficient.h"

#include "lidort_rt.h"
#include "twostream_rt.h"
#include "first_order_rt.h"

namespace FullPhysics {

/****************************************************************//**
*******************************************************************/

class PcaRt : public RadiativeTransferFixedStokesCoefficient {
public:
    PcaRt(const boost::shared_ptr<RtAtmosphere>& Atm,
          const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
          const blitz::Array<double, 1>& Sza, 
          const blitz::Array<double, 1>& Zen, 
          const blitz::Array<double, 1>& Azm,
          int Number_streams, 
          int Number_moments, 
          bool do_solar_sources = true, 
          bool do_thermal_emission = false);

    virtual ~PcaRt() = default;

    virtual int number_stream() const = 0;

    virtual int number_stokes() const { return stokes_coef->stokes_coefficient().cols(); }

    virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
    virtual ArrayAd<double, 2> stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const;

    const boost::shared_ptr<RtAtmosphere>& atmosphere() const { return atm; }
  
    virtual void print(std::ostream& Os, bool Short_form = false) const;

private:

  boost::shared_ptr<RtAtmosphere> atm;

  boost::shared_ptr<LidortRt> lidort_rt;
  boost::shared_ptr<TwostreamRt> twostream_rt;
  boost::shared_ptr<FirstOrderRt> first_order_rt;
};
}
#endif
