#ifndef MAPPING_IMP_BASE_H
#define MAPPING_IMP_BASE_H

#include <string>
#include "array_ad.h"
#include "generic_object.h"

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum
  interface.

  However, almost always you will want to derive from this class
  instead. See PressureImpBase for a more complete discussion of this.

  Individual components of the State Vector may have different
  representations when being retrieved/perturbed rather than being
  used in forward model calculations (e.g. to add a scaling
  factor or log encode the values).

  Derived classes capture those different representations.
*******************************************************************/
class MappingImpBase : public virtual GenericObject {
public:
    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const = 0;

    //-----------------------------------------------------------------------
    /// Calculation of retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const = 0;

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() = 0;

    virtual ~MappingImpBase() {};

};
}

#endif
