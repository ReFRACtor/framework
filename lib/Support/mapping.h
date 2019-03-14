#ifndef MAPPING_H
#define MAPPING_H

#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  Individual SubStateVectorArray components of a State Vector
  may want to store their coefficients in one way, but expose
  them to the retrieval class in another (e.g. to add a scaling
  factor or log encode the values).

  This class captures those mappings from stored coefficients to
  coefficients provided to retrieval processes.

  This class implements 1-to-1 mapping and just provides the
  SubStateVectorArray's coefficients as-is.

  Derived classes can be implemented to extend the mapping concept
  to other mapping types (e.g. scaled, log, etc)
*******************************************************************/
class Mapping {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    Mapping() {};

    //-----------------------------------------------------------------------
    /// Application of mapping (1-to-1 mapping by default)
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1>& apply(ArrayAd<double, 1> const& coeff) const { return coeff; }

    //-----------------------------------------------------------------------
    /// Inversion of mapping (1-to-1 inversion by default)
    //-----------------------------------------------------------------------

    virtual AutoDerivative<double> invert(AutoDerivative<double> coeff_i) const { return coeff_i; }


    virtual ~Mapping() {};
};
}
#endif
