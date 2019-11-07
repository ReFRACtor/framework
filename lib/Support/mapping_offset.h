#ifndef MAPPING_OFFSET_H
#define MAPPING_OFFSET_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a offset for the retrieval view while
  applying that offset to get the forward model view.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingOffset : public Mapping {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingOffset(double Offset, blitz::Array<double, 1> Offsetee)
    : map_name("Offset"), initial_offset(Offset), offsetee(Offsetee) {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff) const {
        int nvar = updated_coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(offsetee + updated_coeff.value()(0)), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const {
        blitz::Array<double, 1> val(1);
        val(0) = initial_offset;
        return ArrayAd<double, 1>(val, true);
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual boost::shared_ptr<Mapping> clone() const
    {
      return boost::shared_ptr<Mapping>(new MappingOffset(initial_offset, offsetee));
    }

    virtual ~MappingOffset() {};

private:
    std::string map_name;
    double initial_offset;
    blitz::Array<double, 1> offsetee; // values being offset

};
}

#endif
