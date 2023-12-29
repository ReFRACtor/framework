#ifndef XSEC_TABLE_IMP_BASE_H
#define XSEC_TABLE_IMP_BASE_H

#include "xsec_table.h"

#include "array_with_unit.h"
#include "linear_interpolate.h"

namespace FullPhysics {

/****************************************************************//**
 Contains common implementation functionality for cross section
 tables.
*******************************************************************/

class XSecTableImpBase: public XSecTable {
public:

  //-----------------------------------------------------------------------
  /// Construct a cross section table instance using a spectral grid
  /// of wavelengths or wavenumbers associated with the cross section values.
  /// The conversion factor is applied after interpolation.
  //-----------------------------------------------------------------------
  
  XSecTableImpBase(const ArrayWithUnit<double, 1>& Spectral_grid,
		   const ArrayWithUnit<double, 2>& XSec_values,
		   double Conversion_factor);

  //-----------------------------------------------------------------------
  /// Get the spectral grid associated with the cross section values
  //-----------------------------------------------------------------------

  virtual const ArrayWithUnit<double, 1> spectral_grid() const
  { return spectral_grid_values; }

  //-----------------------------------------------------------------------
  /// Return the table cross section value at the given spectral point
  /// The value returned is interpolated from the table's grid. A conversion
  /// factor is applied.
  //-----------------------------------------------------------------------

  virtual const DoubleWithUnit cross_section_value(DoubleWithUnit& spectral_point) const;

  //-----------------------------------------------------------------------
  /// Returns additional coefficients such as temperature dependant
  /// coefficients. The value is interpolated from the table's spectral
  /// grid and a conversion factor applied.
  //-----------------------------------------------------------------------

  virtual const ArrayWithUnit<double, 1> cross_section_coefficients(DoubleWithUnit& spectral_point) const;

protected:

  XSecTableImpBase(const ArrayWithUnit<double, 1>& Spectral_grid, const std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > Data_interp, const Unit& Data_units, double Conversion_factor);

  ArrayWithUnit<double, 1> spectral_grid_values;
  
  double conversion_factor;
  Unit xsec_units;

  std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > data_interp;

  XSecTableImpBase() = default;

private:

  void init_interpolation(const ArrayWithUnit<double, 2>& xsec_values);

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(XSecTableImpBase);

#endif
