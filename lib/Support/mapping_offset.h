#ifndef MAPPING_OFFSET_H
#define MAPPING_OFFSET_H

#include <blitz/array.h>

#include "array_ad.h"
#include "mapping_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a offset for the retrieval view while
  applying that offset to get the forward model view.

  For additional information see docs for MappingImpBase class.
*******************************************************************/
class MappingOffset : public MappingImpBase {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingOffset(double Offset, blitz::Array<double, 1> Offsetee)
    : map_name("Offset"), initial_offset(Offset), offsetee(Offsetee) {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& coeff) const {
        int nvar = coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(offsetee + coeff.value()(0)), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_view(ArrayAd<double, 1> const& coeff) const {
        blitz::Array<double, 1> val(1);
        val(0) = initial_offset;
        return ArrayAd<double, 1>(val, true);
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual ~MappingOffset() {};

private:
    std::string map_name;
    double initial_offset;
    blitz::Array<double, 1> offsetee; // values being offset

};
}

#endif
