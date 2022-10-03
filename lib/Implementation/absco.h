#ifndef ABSCO_H
#define ABSCO_H
#include "gas_absorption.h"
#include "array_ad_with_unit.h"
#include "array_with_unit.h"
#include "null_deleter.h"
#include <boost/any.hpp>
#include <vector>

namespace FullPhysics {
class Absco;

/****************************************************************//**
  This is a helper class that calculates the absorption cross section
  for a fixed set of Pressure, Temperature, and Broadener VMR (e.g. H2O
  VMR). It turns out this is a bottle neck, which we can speed up if
  we know ahead of time the values to use.
*******************************************************************/

class AbscoInterpolator : public Printable<AbscoInterpolator> {
public:
  AbscoInterpolator(const boost::shared_ptr<Absco>& A,
		    const ArrayWithUnit<double, 1>& Press, 
		    const ArrayAdWithUnit<double, 1>& Temp,
		    const ArrayAdWithUnit<double, 2>& Broadener_vmr);
  virtual ~AbscoInterpolator() {}
  virtual void print(std::ostream& Os) const {Os << "AbscoInterpolator";}

  blitz::Array<double, 1> absorption_cross_section_noderiv(double wn) const;
  ArrayAd<double, 1> absorption_cross_section_deriv(double wn) const;
private:
  double interpol(double X, bool is_reversed, const std::vector<double>& Xv, 
		  int& i, double& df_dx) const;
  template<class T> blitz::Array<double, 1> 
  absorption_cross_section_noderiv_calc(double wn) const;
  template<class T> ArrayAd<double, 1> 
  absorption_cross_section_deriv_calc(double wn) const;
  inline bool have_no_broadner() const;
  bool have_one_broadner() const
  { return !have_no_broadner() && fb.rows() == 1;}
  bool have_two_broadner() const
  { return !have_no_broadner() && fb.rows() == 2;}
  boost::shared_ptr<Absco> absco;
  blitz::Array<double, 1> p;     // This is in Pascals
  bool p_reversed;		 // True if p goes from high to low
				 // rather than low to high
  ArrayAd<double, 1> t;		// This is in K
  ArrayAd<double, 2> b;		// This is dimensionless
  blitz::Array<double, 2> t_jac; // t.jacobian(), in C order and contiguous
  blitz::Array<double, 3> b_jac; // b.jacobian(), in C order and contiguous
  // Various indexes used in the interpolation
  blitz::Array<int, 1> ip, itp1, itp2;
  blitz::Array<int, 2> ib, ib2;
  blitz::Array<double, 1> dftp1_dt, dftp2_dt, fp, ftp1, ftp2;
  blitz::Array<double, 2> dfb_db, fb;
  mutable ArrayAd<double, 1> res; // We return the same results array
				// each time, so we don't need to keep
				// creating and destroying this.
};

/****************************************************************//**
  This class is used to read the absco tables. 
*******************************************************************/
class Absco: public GasAbsorption {
public:
  virtual ~Absco() {}

//-----------------------------------------------------------------------
/// Scale to apply to underlying ABSCO data to get the 
/// absorption_cross_section. This allows empirical corrections to be
/// applied the ABSCO tables (e.g., O2 scaling). Note that as a user
/// you don't need to apply this correction, it is already applied in
/// absorption_cross_section() and AbscoInterpolator.
//-----------------------------------------------------------------------

  virtual double table_scale(double wn) const = 0;

//-----------------------------------------------------------------------
/// Number of broadener VMR values in absco file. This may be 0, if 
/// we don't have any broadening.
//-----------------------------------------------------------------------

  virtual int number_broadener_vmr(int Broadener_index) const 
  { return broadener_vmr_grid(Broadener_index).rows(); }

//-----------------------------------------------------------------------
/// Number of pressure layers in absco file.
//-----------------------------------------------------------------------

  virtual int number_layer() const 
  {return pressure_grid().rows(); }

//-----------------------------------------------------------------------
/// Number of temperature values in absco file.
//-----------------------------------------------------------------------

  virtual int number_temperature() const 
  {return temperature_grid().cols(); }

//-----------------------------------------------------------------------
/// Return the broadener VMR grid used for this Absco file. This is
/// number_broadener_vmr(Broadener_index) in size, which may be size 0.
///
/// This is dimensionless. 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1>
  broadener_vmr_grid(int Broadener_index) const = 0;

//-----------------------------------------------------------------------
/// Return the pressure grid used for this Absco file. This is
/// number_layer() in size.
///
/// This is in Pascals.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> pressure_grid() const = 0;

//-----------------------------------------------------------------------
/// Return the temperature grid for this Absco file. This is
/// number_layer() x number_temperature() in size. 
///
/// This is in Kelvin.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> temperature_grid() const = 0;


  virtual DoubleWithUnit absorption_cross_section
  (double Wn, const DoubleWithUnit& Press, const DoubleWithUnit& Temp,
    const ArrayWithUnit<double, 1>& Broadener_vmr) const;
  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double wn, 
    const DoubleWithUnit& press, 
    const AutoDerivativeWithUnit<double>& temp,
    const ArrayAdWithUnit<double, 1>& Broadener_vmr) const;

//-----------------------------------------------------------------------
/// Note some of the Absco data files used float, and some use
/// doubles. HDF is happy to convert for us, but because the absco
/// tables are a bit of a bottle neck we actually want to keep track
/// of float vs. double and use this to optimize the read functions.
//-----------------------------------------------------------------------

  virtual bool is_float() const = 0;

//-----------------------------------------------------------------------
/// Read the record for the highest wavenumber that is less than or
/// equal to given wave number. The value returned in number_layer() x
/// number_temperature() x max(number_broadener_vmr(), 1)
///
/// Different absco files have different windows of data. If you don't
/// know the file contains data for wn, you should first check if can
/// read (e.g., use have_data).
///
/// Note that the older binary and 3D Absco files do not have
/// broadener VMR dimension. As a convenience, we pretend like there
/// is a dimension with a single value. This allows us to treat 3D and
/// 4D tables the same.
///
/// \param wn Wave number to read
/// \return Absorption cross section in cm^2 / molecule
//-----------------------------------------------------------------------
  template <class T, int D> blitz::Array<T, D> read(double wn) const;
  friend class AbscoInterpolator;

//-----------------------------------------------------------------------
/// Return either a blitz::Array<double, 3> or blitz::Array<float, 3>, 
/// depending on the type of the underlying data.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 3> read_double(double wn) const = 0;
  virtual blitz::Array<float, 3> read_float(double wn) const = 0;
  virtual blitz::Array<double, 4> read_double_2b(double UNUSED(wn)) const
  { throw Exception("Don't support 2 broadners"); }
  virtual blitz::Array<float, 4> read_float_2b(double UNUSED(wn)) const 
  { throw Exception("Don't support 2 broadners"); }

private:
  double interpol(double X, const std::vector<double>& Xv, 
		  int& i, double& df_dx) const;
  mutable std::vector<double> pgrid;
  mutable std::vector<std::vector<double> > tgrid;
  mutable std::vector<int> tstart_offset;
  mutable std::vector<std::vector<double> > bgrid;
  void fill_pgrid_tgrid_and_bgrid() const;  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

template <> inline blitz::Array<double, 3> Absco::read<double, 3>(double wn) const
{
  return read_double(wn);
}

template <> inline blitz::Array<float, 3> Absco::read<float, 3>(double wn) const
{
  return read_float(wn);
}

template <> inline blitz::Array<double, 4> Absco::read<double, 4>(double wn) const
{
  return read_double_2b(wn);
}

template <> inline blitz::Array<float, 4> Absco::read<float, 4>(double wn) const
{
  return read_float_2b(wn);
}

// We need to define this down here because the definition of Absco
// needs to be defined before we can use it.
inline bool AbscoInterpolator::have_no_broadner() const { return absco->bgrid.size() ==0; }
  
}

FP_EXPORT_KEY(Absco);
#endif

