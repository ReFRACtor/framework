#ifndef MAPPING_SCALE_H
#define MAPPING_SCALE_H

#include <blitz/array.h>

#include "array_ad.h"
#include "mapping_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a scale factor for the retrieval view
  while applying that scale factor to get the forward model view.

  For additional information see docs for MappingImpBase class.
*******************************************************************/
class MappingScale : public MappingImpBase  {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingScale(double Scale, blitz::Array<double, 1> Scalee)
    : map_name("Scale"), initial_scale_factor(Scale), scalee(Scalee) {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const {
        int nvar = coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(scalee * coeff.value()(0)), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const {
        blitz::Array<double, 1> val(1);
        val(0) = initial_scale_factor;
        return ArrayAd<double, 1>(val, true);
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual ~MappingScale() {};

private:
    std::string map_name;
    double initial_scale_factor;
    blitz::Array<double, 1> scalee; // values being scaled by scale factor

};
}

#endif
