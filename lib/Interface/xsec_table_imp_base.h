#ifndef XSEC_TABLE_H
#define XSEC_TABLE_H

#include "xsec_table.h"

namespace FullPhysics {

/****************************************************************//**
 Contains common implementation functionality for cross section
 tables.
*******************************************************************/

class XSecTableFileBase: public XSecTable {

public:

    XSecTableFileBase(const ArrayWithUnit<double, 1>& Spectral_grid, const Array<double, 1>& XSec_values, double Conversion_factor);

    virtual const ArrayWithUnit<double, 1> spectral_grid() const { return spectral_grid_values; }

    virtual const double cross_section_value(DoubleWithUnit& spectral_point) const;

private:

    ArrayWithUnit<double, 1> spectral_grid_values;

    double conversion_factor;

    std::vector<LinearInterpolate<double, double> > data_interp;
}

#endif
