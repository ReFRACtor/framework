#ifndef RAYLEIGH_YOUNG_H
#define RAYLEIGH_YOUNG_H

#include "rayleigh_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the Rayleigh optical depth using the
  Young 1980 method for determining the depolarization factor
  used in the computation of Rayleigh cross section values.

  Young, A.T. (1980), Revised depolarization corrections for atmospheric extinction, 
  Applied Optics, 19 (20), 3427-3428.
*******************************************************************/
class RayleighYoung: public RayleighImpBase {
public:
    RayleighYoung(const boost::shared_ptr<Pressure>& Pres,
                  const std::vector<boost::shared_ptr<Altitude> >& Alt,
                  const boost::shared_ptr<Constant>& C);

    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const;

    virtual boost::shared_ptr<Rayleigh> clone() const;

    virtual void print(std::ostream& Os) const
    {
        Os << "RayleighYoung";
    }

private:
  RayleighYoung() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(RayleighYoung);
#endif
