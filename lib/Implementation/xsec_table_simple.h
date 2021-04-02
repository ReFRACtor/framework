#ifndef XSEC_TABLE_SIMPLE_H
#define XSEC_TABLE_SIMPLE_H

#include "xsec_table_imp_base.h"

namespace FullPhysics {

/****************************************************************//**
 Implements computation of cross section values from a table that
 does not include temperature dependance. Hence there is only a 
 set of directly applicable cross section values in the table.
*******************************************************************/

class XSecTableSimple: public XSecTableImpBase {
public:
    XSecTableSimple(const ArrayWithUnit<double, 1>& Spectral_grid, const blitz::Array<double, 2>& XSec_values, double Conversion_factor);

    virtual ArrayAd<double, 1> optical_depth_each_layer_unweighted(DoubleWithUnit spectral_point, ArrayAd<double, 1> gas_density_levels, ArrayAd<double, 1> temperature_levels) const;

    virtual boost::shared_ptr<XSecTable> clone() const;

private:
    XSecTableSimple(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, double Conversion_factor);

    XSecTableSimple();
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(XSecTableSimple);

#endif
