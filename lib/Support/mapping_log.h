#ifndef MAPPING_LOG_H
#define MAPPING_LOG_H

#include <blitz/array.h>

#include "array_ad.h"
#include "mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements log encoding of coeffs for the retrieval view
  while using the real (linear) values for the forward model view.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingLog : public Mapping {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingLog() : map_name("Log") {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const {
        int nvar = coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(exp(coeff.value())), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const {
        int nvar = coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(log(coeff.value())), nvar, true);
    };

    virtual ~MappingLog() {};

private:
    std::string map_name;

};
}

#endif
