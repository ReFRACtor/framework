#ifndef MAPPING_H
#define MAPPING_H


#include <string>

#include "array_ad.h"
#include "generic_object.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class manages mapping between the State Vector and forward
  model.

  Individual components of the State Vector may have different
  representations (e.g. to add a scaling factor or log encode the
  values) when being retrieved/perturbed vs. the representations
  used in forward model calculations.

  Derived classes capture those different representations.
*******************************************************************/
class Mapping : public virtual GenericObject {
public:
    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const = 0;

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const = 0;

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const = 0;

    virtual boost::shared_ptr<Mapping> clone() const = 0;

    virtual ~Mapping() {};

};
}

#endif
