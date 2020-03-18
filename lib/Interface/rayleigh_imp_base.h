#ifndef RAYLEIGH_IMP_BASE_H
#define RAYLEIGH_IMP_BASE_H

#include "rayleigh.h"

#include <vector>

#include "double_with_unit.h"
#include "constant.h"

#include "pressure.h"
#include "altitude.h"

namespace FullPhysics {

/****************************************************************//**
  This class contains common behavior for Rayleigh optical depth
  implementations. Inheriting classes need only implement the method
  for computing cross section values.
 *******************************************************************/

class RayleighImpBase: public Observer<Pressure>, public Observer<Altitude>, public Rayleigh {

public:

    RayleighImpBase(const boost::shared_ptr<Pressure>& Pres, 
                    const std::vector<boost::shared_ptr<Altitude> >& Alt,
                    const Constant& C);

    virtual void notify_update(const Pressure& UNUSED(P))
    {
        cache_is_stale = true;
    }
    virtual void notify_update(const Altitude& UNUSED(A))
    {
        cache_is_stale = true;
    }

    virtual ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const;

    virtual DoubleWithUnit cross_section(const DoubleWithUnit& W) const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "RayleighImpBase";
    }

private:

    boost::shared_ptr<Pressure> pres;
    std::vector<boost::shared_ptr<Altitude> > alt;

    mutable bool cache_is_stale;

    mutable ArrayAd<double, 2> part_independent_wn;
    void fill_cache() const;

    // Constants. We get this from the Constant class, but stash a copy
    // of them here.
    double molar_weight_dry_air;
};
}
#endif
