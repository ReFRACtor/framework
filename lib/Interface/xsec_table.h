#ifndef XSEC_TABLE_H
#define XSEC_TABLE_H

#include <boost/shared_ptr.hpp>

#include "printable.h"
#include "array_ad_with_unit.h"
#include "double_with_unit.h"

namespace FullPhysics {

/****************************************************************//**
 Represents a cross section table that can compute optical depths
 given its internal representation given density and temperature.
 Note that tables that are not dependant on temperature can safely
 ignore the values.
*******************************************************************/

class XSecTable: public Printable<XSecTable> {
public:

    //-----------------------------------------------------------------------
    /// Compute the unweighted optical depth for each layer at the given 
    /// spectral point given density and temperature on levels.
    ///
    /// Unweighted here means that it is the average of the two levels, 
    /// but has not had a height or pressure difference weighting applied.
    /// 
    /// Not all implementations may have a temperature dependence.
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const = 0;

    //-----------------------------------------------------------------------
    /// Clone the object into a new copy
    //-----------------------------------------------------------------------

    virtual boost::shared_ptr<XSecTable> clone() const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "XSecTable";
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(XSecTable);

#endif
