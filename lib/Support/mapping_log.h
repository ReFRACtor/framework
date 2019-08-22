#ifndef MAPPING_LOG_H
#define MAPPING_LOG_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "mapping_imp_base.h"


namespace FullPhysics {
/****************************************************************//**
  This class implements log encoding of coeffs for the retrieval view
  while using the real (linear) values for the forward model view.

  For additional information see docs for MappingImpBase class.
*******************************************************************/
class MappingLog : public MappingImpBase  {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    MappingLog() : map_name("Log") {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff,
            const boost::shared_ptr<Pressure>& updated_press) const {
        int nvar = updated_coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(exp(updated_coeff.value())), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const {
        int nvar = initial_coeff.number_variable();
        return ArrayAd<double, 1>( blitz::Array<double, 1>(log(initial_coeff.value())), nvar, true);
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual boost::shared_ptr<MappingImpBase> clone() const
    {
      return boost::shared_ptr<MappingImpBase>(new MappingLog());
    }

    virtual ~MappingLog() {};

private:
    std::string map_name;

};
}

#endif
