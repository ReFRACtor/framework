#ifndef MAPPING_LOG_H
#define MAPPING_LOG_H

#include <blitz/array.h>

#include "array_ad.h"
#include "mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements log mapping for SubStateVectorArray
  coefficients.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingLog : public Mapping {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingLog() {};

    //-----------------------------------------------------------------------
    /// Application of log mapping
    /// TODO: Currently, is_const = true. Call number_variable() & diff constr.
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> apply(ArrayAd<double, 1> const& coeff) const {
        return ArrayAd<double, 1>( blitz::Array<double, 1>(exp(coeff.value())), true );
    }

    //-----------------------------------------------------------------------
    /// Application of log mapping
    //-----------------------------------------------------------------------

    virtual const blitz::Array<double, 1> apply(blitz::Array<double, 1> const& coeff) const {
        return blitz::Array<double, 1>(exp(coeff));
    }

    //-----------------------------------------------------------------------
    /// Application of log mapping for a single element
    //-----------------------------------------------------------------------

    virtual AutoDerivative<double> apply_element(AutoDerivative<double> coeff_i) const {
        return exp(coeff_i.value());
    }

    //-----------------------------------------------------------------------
    /// Inversion of log mapping
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> invert(ArrayAd<double, 1> const& coeff) const {
        return ArrayAd<double, 1>( blitz::Array<double, 1>(log(coeff.value())), true);
    }

    //-----------------------------------------------------------------------
    /// Inversion of log mapping
    //-----------------------------------------------------------------------

    virtual const blitz::Array<double, 1> invert(blitz::Array<double, 1> const& coeff) const {
        return blitz::Array<double, 1>(log(coeff));
    }

    //-----------------------------------------------------------------------
    /// Inversion of log mapping for a single element
    //-----------------------------------------------------------------------

    virtual AutoDerivative<double> invert_element(AutoDerivative<double> coeff_i) const {
        return log(coeff_i.value());
    }


    virtual ~MappingLog() {};
};
}
#endif
