#ifndef MAPPING_SCALE_H
#define MAPPING_SCALE_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

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

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff,
            const boost::shared_ptr<Pressure>& updated_press) const {
        int nvar = updated_coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(scalee * updated_coeff.value()(0)), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& updated_coeff) const {
        blitz::Array<double, 1> val(1);
        val(0) = initial_scale_factor;
        return ArrayAd<double, 1>(val, true);
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual boost::shared_ptr<MappingImpBase> clone() const
    {
      return boost::shared_ptr<MappingImpBase>(new MappingScale(initial_scale_factor, scalee));
    }

    virtual ~MappingScale() {};

private:
    std::string map_name;
    double initial_scale_factor;
    blitz::Array<double, 1> scalee; // values being scaled by scale factor

};
}

#endif
