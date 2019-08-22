#ifndef MAPPING_H
#define MAPPING_H

#include "mapping_imp_base.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements 1-to-1 mapping and just provides the
  SubStateVectorArray's coefficients as-is.

  For additional information see docs for MappingImpBase class.
*******************************************************************/
class Mapping : public MappingImpBase {
public:
    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    Mapping() : map_name("OneToOne") {};

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff,
            const boost::shared_ptr<Pressure>& updated_press) const { return updated_coeff; };

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const { return initial_coeff; };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const { return map_name; }

    virtual boost::shared_ptr<MappingImpBase> clone() const
    {
      return boost::shared_ptr<MappingImpBase>(new Mapping());
    }


    virtual ~Mapping() {};

private:
    std::string map_name;
};
}

#endif
