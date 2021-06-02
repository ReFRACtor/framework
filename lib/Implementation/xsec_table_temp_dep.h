#ifndef XSEC_TABLE_TEMP_DEP_H
#define XSEC_TABLE_TEMP_DEP_H

#include "xsec_table_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
 Implements computation of cross section values from a table that
 is dependent on temperature values and therefore uses extra
 coefficient from the cross section table.

 Currently the only known temperature dependent cross section table
 available is for ozone.
*******************************************************************/

class XSecTableTempDep: public XSecTableImpBase {
public:
    XSecTableTempDep(const ArrayWithUnit<double, 1>& Spectral_grid, const ArrayWithUnit<double, 2>& XSec_values, double Conversion_factor);

    virtual ArrayAdWithUnit<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAdWithUnit<double, 1> gas_density_levels, ArrayAdWithUnit<double, 1> temperature_levels) const;

    virtual boost::shared_ptr<XSecTable> clone() const;

private:
    XSecTableTempDep(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, const Unit& Data_units, double Conversion_factor);

    XSecTableTempDep() = default;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(XSecTableTempDep);

#endif
