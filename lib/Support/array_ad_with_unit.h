#ifndef ARRAY_AD_WITH_UNIT_H
#define ARRAY_AD_WITH_UNIT_H
#include "array_ad.h"
#include "auto_derivative_with_unit.h"
#include "unit.h"

namespace FullPhysics {
/****************************************************************//**
  This is a ArrayAd that also has units associated with it.
  This is a simple structure that just keeps these two things
  together. 
*******************************************************************/

template<class T, int D> class ArrayAdWithUnit :
    public Printable<ArrayAdWithUnit<T, D> >
{
public:
  ArrayAdWithUnit() {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const Unit& U)
    : value(V), units(U) {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V, const std::string& U)
    : value(V), units(U) {}
  ArrayAdWithUnit(const ArrayAd<T, D>& V)
    : value(V), units(units::dimensionless) {}
  ArrayAdWithUnit(const ArrayAdWithUnit<T, D>& V)
  {
    reference(V);
  }
  ArrayAdWithUnit& operator=(const ArrayAdWithUnit<T, D>& V)
  {
    reference(V);
    return *this;
  }
  virtual ~ArrayAdWithUnit() = default;
  void print(std::ostream& Os) const 
  { Os << "ArrayAdWithUnit"; }
  ArrayAd<T, D> value;
  Unit units;

//-----------------------------------------------------------------------
/// Convert to the given units.
//-----------------------------------------------------------------------

  inline ArrayAdWithUnit<T, D> convert(const Unit& R) const
  { ArrayAd<T, D> res(value.copy());
    double c = FullPhysics::conversion(units, R);
    res.value() *= c;
    res.jacobian() *= c;
    return ArrayAdWithUnit<T, D>(res, R);
  }
  AutoDerivativeWithUnit<T> operator()(int i1) const
  { return AutoDerivativeWithUnit<T>(value(i1), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3, int i4) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3, i4), units); }
  AutoDerivativeWithUnit<T> operator()(int i1, int i2, int i3, int i4, int i5) const
  { return AutoDerivativeWithUnit<T>(value(i1, i2, i3, i4, i5), units); }
  ArrayAdWithUnit<T, 1> operator()(const blitz::Range& r1) const
  {
    return ArrayAdWithUnit<T, 1>(value(r1), units);
  }
  template<typename T1, typename T2>
  ArrayAdWithUnit<T,blitz::SliceInfo<T,T1,T2>::rank> operator()(T1 r1, T2 r2) const
  {
    return ArrayAdWithUnit<T, blitz::SliceInfo<T,T1,T2>::rank>
      (value(r1, r2), units);
  }
  template<typename T1, typename T2, typename T3>
  ArrayAdWithUnit<T,blitz::SliceInfo<T,T1,T2,T3>::rank> 
  operator()(T1 r1, T2 r2, T3 r3) const
  {
    return ArrayAdWithUnit<T, blitz::SliceInfo<T,T1,T2,T3>::rank>
      (value(r1, r2, r3), units);
  }
  template<typename T1, typename T2, typename T3, typename T4>
  ArrayAdWithUnit<T,blitz::SliceInfo<T,T1,T2,T3,T4>::rank> 
  operator()(T1 r1, T2 r2, T3 r3, T4 r4) const
  {
    return ArrayAdWithUnit<T, blitz::SliceInfo<T,T1,T2,T3,T4>::rank>
      (value(r1, r2, r3, r4), units);
  }
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  ArrayAdWithUnit<T,blitz::SliceInfo<T,T1,T2,T3,T4,T5>::rank> 
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5) const
  {
    return ArrayAdWithUnit<T, blitz::SliceInfo<T,T1,T2,T3,T4,T5>::rank>
      (value(r1, r2, r3, r4, r5), units);
  }
  int rows() const {return value.rows();}
  int cols() const {return value.cols();}
  int depth() const {return value.depth();}
  bool is_constant() const {return value.is_constant();}
  int number_variable() const {return value.number_variable(); }
  void reference(const ArrayAdWithUnit<T, D>& V)
  {
    value.reference(V.value);
    units = V.units;
  }
  inline ArrayAdWithUnit<T, D> copy() const
  { return ArrayAdWithUnit<T, D>(value.copy(), units); }
private:  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

typedef ArrayAdWithUnit<double, 1>  ArrayAdWithUnit_double_1;
typedef ArrayAdWithUnit<double, 2>  ArrayAdWithUnit_double_2;
typedef ArrayAdWithUnit<double, 3>  ArrayAdWithUnit_double_3;
typedef ArrayAdWithUnit<double, 4>  ArrayAdWithUnit_double_4;
}

FP_EXPORT_KEY(ArrayAdWithUnit_double_1);
FP_EXPORT_KEY(ArrayAdWithUnit_double_2);
FP_EXPORT_KEY(ArrayAdWithUnit_double_3);
FP_EXPORT_KEY(ArrayAdWithUnit_double_4);
#endif
