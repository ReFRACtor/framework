#ifndef MAPPING_H
#define MAPPING_H

#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  Individual components of the State Vector may have different
  representations when being retrieved/perturbed rather than being
  used in forward model calculations (e.g. to add a scaling
  factor or log encode the values).

  This class and its subclasses capture those different representations.

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

    Mapping() : map_name("OneToOne") {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const { return coeff; };

    //-----------------------------------------------------------------------
    /// Calculation of retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const { return coeff; };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual ~Mapping() {};

private:
    std::string map_name;
};
}

#endif
