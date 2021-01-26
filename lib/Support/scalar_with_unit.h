#ifndef SCALAR_WITH_UNIT_H
#define SCALAR_WITH_UNIT_H
#include "printable.h"
#include "unit.h"
#include <cmath>

namespace FullPhysics {
  class SpectralDomain;
/****************************************************************//**
  We frequently have a scalar with units associated with it. This is a
  simple structure that just keeps these two things together.
*******************************************************************/
template<class T>
class ScalarWithUnit : public Printable<ScalarWithUnit<T>>,
		       boost::ordered_field_operators<ScalarWithUnit<T>> {
public:
  ScalarWithUnit() {}
  ScalarWithUnit(const ScalarWithUnit<T>& V)
    : value(V.value), units(V.units) {}
  ScalarWithUnit(T V, const Unit& U)
    : value(V), units(U) {}
  ScalarWithUnit(T V, const std::string& U)
    : value(V), units(U) {}
  ScalarWithUnit(T V)
    : value(V), units(units::dimensionless) {}
  T value;
  Unit units;

//-----------------------------------------------------------------------
/// Basic math operators for class.
//-----------------------------------------------------------------------
  inline ScalarWithUnit<T>& operator*=(const ScalarWithUnit<T>& D)
  { value *= D.value; units *= D.units; return *this;}
  inline ScalarWithUnit<T>& operator/=(const ScalarWithUnit<T>& D)
  { value /= D.value; units /= D.units; return *this;}
  inline ScalarWithUnit<T>& operator+=(const ScalarWithUnit<T>& D)
  { value += D.value * FullPhysics::conversion(D.units, units); return *this;}
  inline ScalarWithUnit<T>& operator-=(const ScalarWithUnit<T>& D)
  { value -= D.value * FullPhysics::conversion(D.units, units); return *this;}

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline ScalarWithUnit<T> convert(const Unit& R) const
  { return ScalarWithUnit<T>(value * FullPhysics::conversion(units, R), R); }

//-----------------------------------------------------------------------
/// We often need to handle conversion from wavenumber to/from
/// wavelength. This is either a normal conversion of the units before
/// and after match in the power of length (so cm^-1 to m^-1), or do
/// an inversion. Since we do this often enough, it is worth having a
/// function that handles this logic.
//-----------------------------------------------------------------------
    
  inline ScalarWithUnit<T> convert_wave(const Unit& R) const
  {
    if(units.is_commensurate(R))
      return convert(R);
    else
      return (1.0 / *this).convert(R);
  }

  ScalarWithUnit<T> convert_wave(const Unit& R,
			      const SpectralDomain& Pixel_grid) const;

  void print(std::ostream& Os) const 
  { Os << value << " " << units.name(); }

};

//-----------------------------------------------------------------------
/// \ingroup Miscellaneous
// Compare ScalarWithUnits
///
/// We define <=, >= and > in terms of this operator.
//-----------------------------------------------------------------------
template<class T>
inline bool operator<(const FullPhysics::ScalarWithUnit<T>& A, const FullPhysics::ScalarWithUnit<T>& B)
{ return A.value < B.convert(A.units).value; }

/// We define != in terms of this operator.
//-----------------------------------------------------------------------
template<class T>
inline bool operator==(const FullPhysics::ScalarWithUnit<T>& A, const FullPhysics::ScalarWithUnit<T>& B)
{ return A.value == B.convert(A.units).value; }


}


//-----------------------------------------------------------------------
/// Math functions.
//-----------------------------------------------------------------------
namespace std {
  // Math functions are in std:: namespace.
  template<class T>
  inline FullPhysics::ScalarWithUnit<T> floor(const FullPhysics::ScalarWithUnit<T>& x)
  {
    return FullPhysics::ScalarWithUnit<T>(::floor(x.value), x.units);
  }

  template<class T>
  inline FullPhysics::ScalarWithUnit<T> ceil(const FullPhysics::ScalarWithUnit<T>& x)
  {
    return FullPhysics::ScalarWithUnit<T>(::ceil(x.value), x.units);
  }

  template<class T>
  inline FullPhysics::ScalarWithUnit<T> round(const FullPhysics::ScalarWithUnit<T>& x)
  {
    return FullPhysics::ScalarWithUnit<T>(::round(x.value), x.units);
  }

  template<class T>
  inline FullPhysics::ScalarWithUnit<T> min(const FullPhysics::ScalarWithUnit<T>& x, const FullPhysics::ScalarWithUnit<T>& y)
  {
    return FullPhysics::ScalarWithUnit<T>(std::min(x.value, y.value * FullPhysics::conversion(y.units, x.units)), x.units);
  }

  template<class T>
  inline FullPhysics::ScalarWithUnit<T> max(const FullPhysics::ScalarWithUnit<T>& x, const FullPhysics::ScalarWithUnit<T>& y)
  {
    return FullPhysics::ScalarWithUnit<T>(std::max(x.value, y.value * FullPhysics::conversion(y.units, x.units)), x.units);
  }

}

#endif
