#ifndef FLOAT_WITH_UNIT_H
#define FLOAT_WITH_UNIT_H
#include "printable.h"
#include "unit.h"
#include <cmath>

namespace FullPhysics {
  class SpectralDomain;
/* TODO: Move common behavior between DoubleWithUnit/FloatWithUnit up to templated VariableWithUnit class */
/****************************************************************//**
  In OSS there are floats with units associated with them. This is a
  an exact copy of DoubleWithUnit for floats.
*******************************************************************/
class FloatWithUnit : public Printable<FloatWithUnit>,
		       boost::ordered_field_operators<FloatWithUnit> {
public:
  FloatWithUnit() {}
  FloatWithUnit(float V, const Unit& U)
    : value(V), units(U) {}
  FloatWithUnit(float V, const std::string& U)
    : value(V), units(U) {}
  FloatWithUnit(float V)
    : value(V), units(units::dimensionless) {}
  float value;
  Unit units;

//-----------------------------------------------------------------------
/// Basic math operators for class.
//-----------------------------------------------------------------------
  inline  FloatWithUnit& operator*=(const FloatWithUnit& D)
  { value *= D.value; units *= D.units; return *this;}
  inline FloatWithUnit& operator/=(const FloatWithUnit& D)
  { value /= D.value; units /= D.units; return *this;}
  inline FloatWithUnit& operator+=(const FloatWithUnit& D)
  { value += D.value * FullPhysics::conversion(D.units, units); return *this;}
  inline FloatWithUnit& operator-=(const FloatWithUnit& D)
  { value -= D.value * FullPhysics::conversion(D.units, units); return *this;}

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline FloatWithUnit convert(const Unit& R) const
  { return FloatWithUnit(value * FullPhysics::conversion(units, R), R); }

//-----------------------------------------------------------------------
/// We often need to handle conversion from wavenumber to/from
/// wavelength. This is either a normal conversion of the units before
/// and after match in the power of length (so cm^-1 to m^-1), or do
/// an inversion. Since we do this often enough, it is worth having a
/// function that handles this logic.
//-----------------------------------------------------------------------
    
  inline FloatWithUnit convert_wave(const Unit& R) const
  {
    if(units.is_commensurate(R))
      return convert(R);
    else
      return (1.0 / *this).convert(R);
  }

  FloatWithUnit convert_wave(const Unit& R,
			      const SpectralDomain& Pixel_grid) const;

  void print(std::ostream& Os) const 
  { Os << value << " " << units.name(); }

};

//-----------------------------------------------------------------------
/// \ingroup Miscellaneous
// Compare FloatWithUnits
///
/// We define <=, >= and > in terms of this operator.
//-----------------------------------------------------------------------
inline bool operator<(const FullPhysics::FloatWithUnit& A, const FullPhysics::FloatWithUnit& B)
{ return A.value < B.convert(A.units).value; }

/// We define != in terms of this operator.
//-----------------------------------------------------------------------
inline bool operator==(const FullPhysics::FloatWithUnit& A, const FullPhysics::FloatWithUnit& B)
{ return A.value == B.convert(A.units).value; }


}


//-----------------------------------------------------------------------
/// Math functions.
//-----------------------------------------------------------------------
namespace std {
  // Math functions are in std:: namespace.
  inline FullPhysics::FloatWithUnit floor(const FullPhysics::FloatWithUnit& x)
  {
    return FullPhysics::FloatWithUnit(::floor(x.value), x.units);
  }

  inline FullPhysics::FloatWithUnit ceil(const FullPhysics::FloatWithUnit& x)
  {
    return FullPhysics::FloatWithUnit(::ceil(x.value), x.units);
  }

  inline FullPhysics::FloatWithUnit round(const FullPhysics::FloatWithUnit& x)
  {
    return FullPhysics::FloatWithUnit(::round(x.value), x.units);
  }

  inline FullPhysics::FloatWithUnit min(const FullPhysics::FloatWithUnit& x, const FullPhysics::FloatWithUnit& y)
  {
    return FullPhysics::FloatWithUnit(std::min(x.value, static_cast<float>(y.value * FullPhysics::conversion(y.units, x.units))), x.units);
  }

  inline FullPhysics::FloatWithUnit max(const FullPhysics::FloatWithUnit& x, const FullPhysics::FloatWithUnit& y)
  {
    return FullPhysics::FloatWithUnit(std::max(x.value, static_cast<float>(y.value * FullPhysics::conversion(y.units, x.units))), x.units);
  }
}

#endif
