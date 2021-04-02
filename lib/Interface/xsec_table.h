#ifndef XSEC_TABLE_H
#define XSEC_TABLE_H

#include <boost/shared_ptr.hpp>

#include "printable.h"
#include "array_ad.h"
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

    virtual ArrayAd<double, 2> optical_depth_each_layer(DoubleWithUnit grid_point, ArrayAd<double, 1> gas_density_levels, ArrayAd<double, 1> temperature_levels) const = 0;

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
