#include "absco.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include <algorithm>
#include <functional>
#include <cmath>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Absco, GasAbsorption)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Find interpolation factor in sorted array.
//-----------------------------------------------------------------------

double Absco::interpol(double X, const std::vector<double>& Xv, 
		       int& i, double& df_dx) const
{
  // Find index number of largest Xv that is < X.
  i = (int) (std::lower_bound(Xv.begin(), Xv.end(), X) - Xv.begin()) - 1;
  if(i < 0)
    i = 0;
  if(i > (int) Xv.size() - 2)
    i = (int) Xv.size() - 2;
  df_dx = 1 / (Xv[i + 1] - Xv[i]);
  return (X - Xv[i]) * df_dx ;
}

//-----------------------------------------------------------------------
/// Find interpolation factor in sorted array.
//-----------------------------------------------------------------------

double AbscoInterpolator::interpol(double X, bool is_reversed,
				   const std::vector<double>& Xv, 
				   int& i, double& df_dx) const
{
  // Find index number of largest Xv that is < X.
  if(is_reversed)
    i = (int) (std::lower_bound(Xv.begin(), Xv.end(), X,
				std::greater<double>()) - Xv.begin()) - 1;
  else
    i = (int) (std::lower_bound(Xv.begin(), Xv.end(), X) - Xv.begin()) - 1;
  if(i < 0)
    i = 0;
  if(i > (int) Xv.size() - 2)
    i = (int) Xv.size() - 2;
  df_dx = 1 / (Xv[i + 1] - Xv[i]);
  return (X - Xv[i]) * df_dx ;
}

//-----------------------------------------------------------------------
/// Fill pgrid and tgrid if needed. This is then used by
/// AbscoInterpolator, but we store it in Absco because we only need
/// to calculate this once, regardless of how many AbscoInterpolator
/// we create.
//-----------------------------------------------------------------------

void Absco::fill_pgrid_tgrid_and_bgrid() const
{
  // Store grid in a form that is easier to search.
  if(pgrid.size() == 0) {
    Array<double, 1> pgridg(pressure_grid());
    Array<double, 2> tg(temperature_grid());
    pgrid.resize(pgridg.rows());
    tgrid.resize(tg.rows());
    tstart_offset.resize(tg.rows());
    bgrid.resize(number_broadener());
    for(int i = 0; i < pgridg.rows(); ++i)
      pgrid[i] = pgridg(i);
    for(int i = 0; i < number_broadener(); ++i) {
      Array<double, 1> bgridg(broadener_vmr_grid(i));
      bgrid[i].resize(bgridg.rows());
      for(int j = 0; j < bgridg.rows(); ++j)
	bgrid[i][j] = bgridg(j);
    }
    // Some temperature grids use nan to indicate that we don't have a
    // value, i.e., the temperature is a jagged array. This is the
    // case for AbscoAer. We just chop off the missing data, but we
    // need to keep a count of what we chop off so we know how to look
    // into the cross section data.
    for(int i = 0; i < tg.rows(); ++i) {
      tgrid[i].clear();
      bool first_good_data = true;
      for(int j = 0; j < tg.cols(); ++j) 
	if(isfinite(tg(i,j))) {
	  if(first_good_data) {
	    first_good_data = false;
	    tstart_offset[i] = j;
	  }
	  tgrid[i].push_back(tg(i,j));
	}
    }
  }
}

// See base class for description.
DoubleWithUnit Absco::absorption_cross_section
(double Wn, const DoubleWithUnit& Press, const DoubleWithUnit& Temp,
 const ArrayWithUnit<double, 1>& Broadener_vmr) const
{
  blitz::Range ra = blitz::Range::all();
  ArrayWithUnit<double, 1> p;
  ArrayAdWithUnit<double, 1> t;
  ArrayAdWithUnit<double, 2> b;
  p.units = Press.units;
  p.value.resize(1);
  p.value(0) = Press.value;
  t.units = Temp.units;
  t.value.resize(1, 0);
  t.value(0) = Temp.value;
  b.units = Broadener_vmr.units;
  b.value.resize(1, Broadener_vmr.rows(), 0);
  b.value(0, ra) = Broadener_vmr.value;
  boost::shared_ptr<Absco> self_ptr = to_ptr(const_cast<Absco&>(*this));
  AbscoInterpolator aint(self_ptr, p, t, b);
  return DoubleWithUnit(aint.absorption_cross_section_noderiv(Wn)(0),
			"cm^2");
}

// See base class for description.
AutoDerivativeWithUnit<double> Absco::absorption_cross_section
(double Wn, const DoubleWithUnit& Press, 
 const AutoDerivativeWithUnit<double>& Temp,
 const ArrayAdWithUnit<double, 1>& Broadener_vmr) const
{
  blitz::Range ra = blitz::Range::all();
  ArrayWithUnit<double, 1> p;
  ArrayAdWithUnit<double, 1> t;
  ArrayAdWithUnit<double, 2> b;
  p.units = Press.units;
  p.value.resize(1);
  p.value(0) = Press.value;
  t.units = Temp.units;
  t.value.resize(1, Temp.value.number_variable());
  t.value(0) = Temp.value;
  b.units = Broadener_vmr.units;
  b.value.resize(1, Broadener_vmr.rows(), Broadener_vmr.number_variable());
  b.value(0, ra) = Broadener_vmr.value;
  boost::shared_ptr<Absco> self_ptr = to_ptr(const_cast<Absco&>(*this));
  AbscoInterpolator aint(self_ptr, p, t, b);
  return AutoDerivativeWithUnit<double>
    (aint.absorption_cross_section_deriv(Wn)(0), "cm^2");
}

//-----------------------------------------------------------------------
/// Set up a AbscoInterpolator for the given set up pressure,
/// temperature, and broadner VMR.
///
/// \param A Absco data to use
/// \param Press Set of pressures.
/// \param Temp Set of Temperatures. Should be same size as
///    press. 
/// \param Broadener_vmr Set of Broadner VMR (e.g., H2O VMR). Not all
///   tables will make use of this information.
//-----------------------------------------------------------------------
AbscoInterpolator::AbscoInterpolator
(const boost::shared_ptr<Absco>& A,
 const ArrayWithUnit<double, 1>& Press, 
 const ArrayAdWithUnit<double, 1>& Temp,
 const ArrayAdWithUnit<double, 2>& Broadener_vmr)
  : absco(A),
    ip(Press.rows()), itp1(Press.rows()), itp2(Press.rows()),
    ib(Broadener_vmr.rows(), Press.rows()),
    ib2(Broadener_vmr.rows(), Press.rows()), dftp1_dt(Press.rows()), 
    dftp2_dt(Press.rows()), dfb_db(Broadener_vmr.rows(), Press.rows()),
    fp(Press.rows()), 
    fb(Broadener_vmr.rows(), Press.rows()),
    ftp1(Press.rows()),
    ftp2(Press.rows()),
    res(Press.rows(), std::max(Temp.value.number_variable(), 
			       Broadener_vmr.value.number_variable()))
{
  if(Press.rows() != Temp.rows() ||
     Press.rows() != Broadener_vmr.cols())
    throw Exception("Pressure, Temperature, and Broadener_vmr arrays needs to be the same size");
  p.reference(Press.convert(units::Pa).value);
  t.reference(Temp.convert(units::K).value);
  b.reference(Broadener_vmr.convert(units::dimensionless).value);
  t_jac.reference(to_c_order(t.jacobian()));
  b_jac.reference(to_c_order(b.jacobian()));

  absco->fill_pgrid_tgrid_and_bgrid();
  p_reversed = (absco->pgrid[0] > absco->pgrid[1]);
  for(int i = 0; i < p.rows(); ++i) {
    double unused;
    fp(i) = interpol(p(i), p_reversed, absco->pgrid, ip(i), unused);
    // We might not actually need to interpolate over broadener
    if(absco->bgrid.size() ==0) {
      ib(0, i) = 0;
      ib2(0, i) = 0;
      dfb_db(0, i) = 0;
      fb(0, i) = 1;
    } else {
      for(int j = 0; j < fb.rows(); ++j) {
	fb(j, i) = interpol(b.value()(j,i), false, absco->bgrid[j], ib(j, i),
			    dfb_db(j, i));
	ib2(j, i) = ib(j, i) + 1;
      }
    }
    ftp1(i) = interpol(t.value()(i), false, absco->tgrid[ip(i)], itp1(i),
		       dftp1_dt(i));
    itp1(i) += absco->tstart_offset[ip(i)];
    ftp2(i) = interpol(t.value()(i), false, absco->tgrid[ip(i) + 1], itp2(i), 
		       dftp2_dt(i));
    itp2(i) += absco->tstart_offset[ip(i)+1];
  }
}

template<class T> Array<double, 1> 
AbscoInterpolator::absorption_cross_section_noderiv_calc(double wn) const
{
  if(fb.rows() == 2) {
    std::cerr << "Hi there\n"
	      << fb << "\n"
	      << ib << "\n";
  }
  Array<T, 3> a(absco->read<T, 3>(wn));
  for(int i = 0; i < res.rows(); ++i) {
    double t11 = a(ip(i), itp1(i), ib(0, i)) * (1 - ftp1(i)) + 
      a(ip(i), itp1(i) + 1, ib(0, i)) * ftp1(i);
    double t12 = a(ip(i) + 1, itp2(i), ib(0,i)) * (1 - ftp2(i)) + 
      a(ip(i) + 1, itp2(i) + 1, ib(0,i)) * ftp2(i);
    double t1 = t11 * (1 - fp(i)) + t12 * fp(i);
    double t21 = a(ip(i), itp1(i), ib2(0,i)) * (1 - ftp1(i)) + 
      a(ip(i), itp1(i) + 1, ib2(0,i)) * ftp1(i);
    double t22 = a(ip(i) + 1, itp2(i), ib2(0,i)) * (1 - ftp2(i)) + 
      a(ip(i) + 1, itp2(i) + 1, ib2(0,i)) * ftp2(i);
    double t2 = t21 * (1 - fp(i)) + t22 * fp(i);
    res.value()(i) = t1 * (1 - fb(0,i)) + t2 * fb(0,i);
  }
  // Apply scale
  res.value() *= absco->table_scale(wn);
  return res.value();
}

template  Array<double, 1> AbscoInterpolator::absorption_cross_section_noderiv_calc<double>(double wn) const;
template  Array<double, 1> AbscoInterpolator::absorption_cross_section_noderiv_calc<float>(double wn) const;

//-----------------------------------------------------------------------
/// Return absorption cross section, without derivatives.
/// \param wn wave number
/// \return Absorption cross section in cm^2 / molecule
//-----------------------------------------------------------------------

Array<double, 1> AbscoInterpolator::absorption_cross_section_noderiv
(double wn) const
{
  if(absco->is_float())
    return absorption_cross_section_noderiv_calc<float>(wn);
  else
    return absorption_cross_section_noderiv_calc<double>(wn);
}

template<class T> ArrayAd<double, 1> 
AbscoInterpolator::absorption_cross_section_deriv_calc(double wn) const
{
  Range ra(Range::all());
  Array<T, 3> a(absco->read<T, 3>(wn));
  // Turns out that jacobian calculation is faster if we use pointers. 
  // This is *not* true for things like itp1(i), the speed is the same 
  // (for gcc 4.6, -O2), so we leave these with the clearer expression.
  double *restrict res_jac_p = res.jacobian().data();
  const double* restrict t_jac_p = t_jac.data();
  const double* restrict b_jac_p = b_jac.data();
  for(int i = 0; i < res.rows(); ++i) {
    double t11 = a(ip(i), itp1(i), ib(0,i)) * (1 - ftp1(i)) + 
      a(ip(i), itp1(i) + 1, ib(0,i)) * ftp1(i);
    double dt11_dt = 
      (a(ip(i), itp1(i) + 1, ib(0,i)) - a(ip(i), itp1(i), ib(0,i))) * dftp1_dt(i);
    double t12 = a(ip(i) + 1, itp2(i), ib(0,i)) * (1 - ftp2(i)) + 
      a(ip(i) + 1, itp2(i) + 1, ib(0,i)) * ftp2(i);
    double dt12_dt =
      (a(ip(i) + 1, itp2(i) + 1, ib(0,i)) - a(ip(i) + 1, itp2(i), ib(0,i))) * 
      dftp2_dt(i);
    double t1 = t11 * (1 - fp(i)) + t12 * fp(i);
    double dt1_dt =
      dt11_dt * (1 - fp(i)) + dt12_dt * fp(i);

    double t21 = a(ip(i), itp1(i), ib2(0,i)) * (1 - ftp1(i)) + 
      a(ip(i), itp1(i) + 1, ib2(0,i)) * ftp1(i);
    double dt21_dt = 
      (a(ip(i), itp1(i) + 1, ib2(0,i)) - a(ip(i), itp1(i), ib2(0,i))) * dftp1_dt(i);
    double t22 = a(ip(i) + 1, itp2(i), ib2(0,i)) * (1 - ftp2(i)) + 
      a(ip(i) + 1, itp2(i) + 1, ib2(0,i)) * ftp2(i);
    double dt22_dt =
      (a(ip(i) + 1, itp2(i) + 1, ib2(0,i)) - a(ip(i) + 1, itp2(i), ib2(0,i))) * 
      dftp2_dt(i);
    double t2 = t21 * (1 - fp(i)) + t22 * fp(i);
    double dt2_dt =
      dt21_dt * (1 - fp(i)) + dt22_dt * fp(i);
    res.value()(i) = t1 * (1 - fb(0,i)) + t2 * fb(0,i);
    double dr_db =
      (t2 - t1) * dfb_db(0, i);
    double dr_dt =
      dt1_dt * (1 - fb(0,i)) + dt2_dt * fb(0,i);
    if(!t.is_constant() && !b.is_constant())
      for(int k = 0; k < res.jacobian().cols(); ++k)
	*res_jac_p++ = dr_dt * *t_jac_p++ + dr_db * *b_jac_p++;
    else if(!t.is_constant())
      for(int k = 0; k < res.jacobian().cols(); ++k)
	*res_jac_p++ = dr_dt * *t_jac_p++;
    else if(!b.is_constant())
      for(int k = 0; k < res.jacobian().cols(); ++k)
	*res_jac_p++ = dr_db * *b_jac_p++;
  }
  // Apply scale
  double s = absco->table_scale(wn);
  res.value() *= s;
  res.jacobian() *= s;
  return res;
}

template  ArrayAd<double, 1> AbscoInterpolator::absorption_cross_section_deriv_calc<double>(double wn) const;
template  ArrayAd<double, 1> AbscoInterpolator::absorption_cross_section_deriv_calc<float>(double wn) const;

//-----------------------------------------------------------------------
/// Return absorption cross section, with derivatives.
/// \param wn wave number
/// \return Absorption cross section in cm^2 / molecule
//-----------------------------------------------------------------------

ArrayAd<double, 1> AbscoInterpolator::absorption_cross_section_deriv
(double wn) const
{
  if(absco->is_float())
    return absorption_cross_section_deriv_calc<float>(wn);
  else
    return absorption_cross_section_deriv_calc<double>(wn);
}
