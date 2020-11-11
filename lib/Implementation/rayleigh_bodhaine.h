#ifndef RAYLEIGH_BODHAINE_H
#define RAYLEIGH_BODHAINE_H

#include "rayleigh_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
  USES Bodhaine et al. (1999) parameterization of rayleigh scattering
  Bodhaine et al., 1999, Journ. Atm. Oc. Tech., 16(11), 1854-1861.

  The full formula is:

                  24 PI^3    (n^2-1)^2   (6+3dpf)
  cross-section = ------- *  ---------- * --------
                 wl^4 Ns^2   (n^2+2)^2   (6-7dpf)

 where n is the index of refraction of standard air at some P,T, and Ns
 is the number of molecules per cubic centimeter AT THE SAME P,T. (the
 product of those terms turns out to be independent of P & T, according
 to Lorenz-Lorentz theory); dpf is the depolarization factor, and wl is
 the wavelength in centimeters.

 The result is then the Rayleigh cross section in cm^2 per molecule of
 air.
*******************************************************************/
class RayleighBodhaine: public RayleighImpBase {
public:
    RayleighBodhaine(const boost::shared_ptr<Pressure>& Pres,
                     const std::vector<boost::shared_ptr<Altitude> >& Alt,
                     const boost::shared_ptr<Constant>& C);

    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const;

    virtual boost::shared_ptr<Rayleigh> clone() const;

    virtual void print(std::ostream& Os) const
    {
        Os << "RayleighBodhaine";
    }

private:
  RayleighBodhaine() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(RayleighBodhaine);
#endif
