#include "aerosol_optical.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolOptical::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aerosol)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverAerosolExtinction)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverAerosolProperty)
    & FP_NVP(aext) & FP_NVP(aprop) & FP_NVP(press) & FP_NVP(rh)
    & FP_NVP_(reference_wn);
}

FP_IMPLEMENT(AerosolOptical);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

typedef const boost::shared_ptr<AerosolExtinction>& (AerosolOptical::*a1)(int) const;

typedef const boost::shared_ptr<AerosolProperty>& (AerosolOptical::*a2)(int) const;

REGISTER_LUA_DERIVED_CLASS(AerosolOptical, Aerosol)
.def(luabind::constructor<const std::vector<boost::shared_ptr<AerosolExtinction> >&,
     const std::vector<boost::shared_ptr<AerosolProperty> >&,
     const boost::shared_ptr<Pressure>&,
     const boost::shared_ptr<RelativeHumidity>&>())
.def("number_particle", &AerosolOptical::number_particle)
.def("aerosol_extinction", ((a1) &AerosolOptical::aerosol_extinction))
.def("aerosol_property", ((a2) &AerosolOptical::aerosol_property))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an aerosol.
/// \param Aext Aerosol extinction for each aerosol.
/// \param Aerosol_prop Aerosol properties for each aerosol.
/// \param Press The Pressure object that gives the pressure grid.
/// \param Rh The RelativeHumidity object that gives the relative humidity.
/// \param Reference_wn The wavenumber that Aext is given for. This
///    is optional, the default value matches the reference band given
///    in the ATB.
//-----------------------------------------------------------------------

AerosolOptical::AerosolOptical
(const std::vector<boost::shared_ptr<AerosolExtinction> >& Aext,
 const std::vector<boost::shared_ptr<AerosolProperty> >& Aerosol_prop,
 const boost::shared_ptr<Pressure>& Press,
 const boost::shared_ptr<RelativeHumidity>& Rh,
 double Reference_wn)
: aext(Aext),
  aprop(Aerosol_prop),
  press(Press),
  rh(Rh),
  reference_wn_(Reference_wn),
  cache_is_stale(true),
  nlay(-1),
  nvar(-1)
{
  if((int) aprop.size() != number_particle())
    throw Exception("aprop needs to be size of number_particle()");
  for(int i = 0; i < number_particle(); ++i) {
    aprop[i]->add_observer(*this);
    aext[i]->add_observer(*this);
  }
  press->add_observer(*this);
}

//-----------------------------------------------------------------------
/// Aerosol names, plus the string "total" as the 1st entry. This
/// matches what is returned by
/// aerosol_optical_depth_each_particle_and_total(). This may seem a
/// bit odd, but this is what is expected as output to the HDF file.
///
/// Note, this is returned as a blitz::Array rather
/// than the more natural std::vector because this is what is needed
/// to write this out to HDF.
//-----------------------------------------------------------------------

blitz::Array<std::string, 1> AerosolOptical::aerosol_name_arr() const
{
  blitz::Array<std::string, 1> res((int) aext.size());
  for(int i = 0; i < res.rows(); i++)
    res(i) = aext[i]->aerosol_name();
  return res;
}

//-----------------------------------------------------------------------
/// We cache the part of the extinction_optical_depth_each_layer calculation that
/// is independent of wn.
//-----------------------------------------------------------------------

void AerosolOptical::fill_cache() const
{
  if(!cache_is_stale)
    return;
  // Initial time through, set nvar if we haven't already set this the
  // state vector update (e.g., we call this before doing a state
  // vector update)
  if(nvar < 0) {
    for(int j = 0; j < number_particle(); ++j) {
      nvar = std::max(nvar, aprop[j]->extinction_coefficient_each_layer(reference_wn_).number_variable());
      for(int i = 0; i < press->number_layer(); ++i)
      nvar = std::max(nvar,
                  aext[j]->extinction_for_layer(i).number_variable());
    }
    for(int i = 0; i < press->number_layer(); ++i)
      nvar = std::max(nvar, press->pressure_grid()(i).value.number_variable());
  }
  od_ind_wn.resize(press->number_layer(), number_particle(), nvar);
  for(int i = 0; i < od_ind_wn.rows(); ++i) {
    AutoDerivativeWithUnit<double> delta_press = 
      (press->pressure_grid()(i + 1) - press->pressure_grid()(i));
    AutoDerivative<double> dp = delta_press.convert(units::Pa).value;

    // Resize number of variables in case surface pressure is not retrieved
    if (dp.number_variable() == 0) {
      dp.gradient().resize(nvar);
      dp.gradient() = 0.0;
    }

    for(int j = 0; j < number_particle(); ++j) {
      /// We scale the extinction coefficient return by aprop at the 
      /// reference wave number, so that aext of 1 means the extinction
      /// coefficient for a particle is 1.
      od_ind_wn(i, j) = 1.0 / aprop[j]->extinction_coefficient_each_layer(reference_wn_)(i) *
        dp * aext[j]->extinction_for_layer(i);
    }
  }
  /*
  std::cout << "# Optical depths for each layer value: " << std::endl
          << od_ind_wn.value() << std::endl
          << "# Optical depths for each layer jacobian: " << std::endl
          << od_ind_wn.jacobian() << std::endl;
  */
  nlay = press->number_layer();
  cache_is_stale = false;
}

//-----------------------------------------------------------------------
/// This gives the extinction optical depth for each layer, for the given wave
/// number. Note this only includes the aerosol portion of this,
/// Atmosphere class combines this with Absorbers and rayleigh
/// scattering.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer() x number_particle().
//-----------------------------------------------------------------------

ArrayAd<double, 2> 
AerosolOptical::extinction_optical_depth_each_layer(double wn) const
{
  Range ra(Range::all());
  firstIndex i1; secondIndex i2; thirdIndex i3;
  FunctionTimer ft(timer.function_timer());
  fill_cache();
  ArrayAd<double, 2> res(od_ind_wn.copy());
  for(int i = 0; i < number_particle(); ++i) {
    ArrayAd<double, 1> t = aprop[i]->extinction_coefficient_each_layer(wn);
    if(res.is_constant() && !t.is_constant()) {
      res.resize_number_variable(t.number_variable());
      res.jacobian() = 0;
    }
    Array<double, 1> v(res.value()(ra, i));
    Array<double, 2> jac(res.jacobian()(ra, i, ra));
    v *= t.value();
    if(t.is_constant())
      jac = t.value()(i1) * jac(i1, i2);
    else
      jac = t.value()(i1) * jac(i1,i2) + v(i1) * t.jacobian()(i1, i2);
  }
  return res;
}

//-----------------------------------------------------------------------
/// This gives the aerosol scattering extinction each layer, for the given wave
/// number. Note this only includes the
/// aerosol portion of this, 
/// Atmosphere class combines this with Rayleigh scattering.
///
/// The equation here is:
///     tau_sca = tau_ind  * k_sca,wn
/// Where tau_ind is the independent portion of the optical depth and k_sca
/// is the scattering cofficient at the requested wave number.
///
/// Equivalently this could be written as:
///     tau_sca = tau_ext * k_sca,wn / k_ext,wn
///
/// Which means that the scattering optical depth is also equal to the
/// extinction optical depth times the aerosol single scattering albedo 
/// calculated as k_sca,wn / k_ext,wn.
///
/// The above relationship is important when dealing with jacobians related
/// to aerosol extinction.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer()
//-----------------------------------------------------------------------

ArrayAd<double, 2> 
AerosolOptical::scattering_optical_depth_each_layer(double wn) const
{
  Range ra(Range::all());
  firstIndex i1; secondIndex i2; thirdIndex i3;
  FunctionTimer ft(timer.function_timer());
  fill_cache();
  ArrayAd<double, 2> res(od_ind_wn.copy());
  for(int i = 0; i < number_particle(); ++i) {
    ArrayAd<double, 1> t = aprop[i]->scattering_coefficient_each_layer(wn);
    if(res.is_constant() && !t.is_constant()) {
      res.resize_number_variable(t.number_variable());
      res.jacobian() = 0;
    }
    Array<double, 1> v(res.value()(ra, i));
    Array<double, 2> jac(res.jacobian()(ra, i, ra));
    v *= t.value();
    if(t.is_constant())
      jac = t.value()(i1) * jac(i1, i2);
    else
      jac = t.value()(i1) * jac(i1,i2) + v(i1) * t.jacobian()(i1, i2);
  }
  return res;
}

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol for a single particle. This is
/// number_layer x number_moments x number_scatter
/// \param wn The wave number.
/// \param pindex The particle index.
/// \return Scattering moments for each layer. This is 
///         number_moment + 1 x number_layer() x number scattering
///         matrix elements
//-----------------------------------------------------------------------

ArrayAd<double, 3> AerosolOptical::pf_mom(double wn, int pindex, int nummom, int numscat) const
{
  FunctionTimer ft(timer.function_timer());
  range_check(pindex, 0, number_particle());
  return aprop[pindex]->phase_function_moment_each_layer(wn, nummom, numscat);
}

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol.
/// \param wn The wave number.
/// \param frac_aer This is number_active_layer() x number_particle()
/// \return Scattering moments for each layer. This is 
///         number_moment + 1 x number_layer() x number scattering
///         matrix elements
//-----------------------------------------------------------------------

blitz::Array<double, 3> AerosolOptical::pf_mom(double wn, 
            const blitz::Array<double, 2>& frac_aer) const
{
  FunctionTimer ft(timer.function_timer());
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra = Range::all();
  if(frac_aer.rows() != nlay ||
     frac_aer.cols() != number_particle()) {
    std::stringstream err_msg;
    err_msg << "frac_aer needs to be number_active_layer = "
          << nlay
          << " x number_particle = "
          << number_particle() << ", but is currently "
          << frac_aer.rows() << " x " << frac_aer.cols();
    throw Exception(err_msg.str());
  }
  std::vector<Array<double, 3> > pf;
  int s1 = 0;
  int s2 = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    pf.push_back(pf_mom(wn, j).value());
    s1 = std::max(s1, pf[j].rows());
    s2 = std::max(s2, pf[j].depth());
  }
  blitz::Array<double, 3> res(s1, nlay, s2);
  res = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    Range r1(0, pf[j].rows() - 1);
    Range r2(0, pf[j].depth() - 1);
    res(r1,ra,r2) += frac_aer(ra,j)(i2) * pf[j](i1, i2, i3);
  }
  return res;
}

ArrayAd<double, 3> AerosolOptical::pf_mom(double wn, 
                           const ArrayAd<double, 2>& frac_aer,
                           int nummom, int numscat) const
{
  FunctionTimer ft(timer.function_timer());
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  Range ra = Range::all();
  if(frac_aer.rows() != nlay ||
     frac_aer.cols() != number_particle()) {
    std::stringstream err_msg;
    err_msg << "frac_aer needs to be number_active_layer = "
          << nlay
          << " x number_particle = "
          << number_particle() << ", but is currently "
          << frac_aer.rows() << " x " << frac_aer.cols();
    throw Exception(err_msg.str());
  }

  std::vector<ArrayAd<double, 3> > pf;
  int s1 = 0;
  int s2 = 0;
  int nvar = frac_aer.number_variable();
  for(int j = 0; j < frac_aer.cols(); ++j) {
    pf.push_back(aprop[j]->phase_function_moment_each_layer(wn, nummom, numscat));
    s1 = std::max(s1, pf[j].rows());
    s2 = std::max(s2, pf[j].depth());
    nvar = std::max(nvar, pf[j].number_variable());
    if(!pf[j].is_constant() && !frac_aer.is_constant() &&
       pf[j].number_variable() != frac_aer.number_variable())
      throw Exception("We don't currently have the code working correctly for combining intermediates and state vector derivatives. We'll need to think through how to do this.");
  }
  ArrayAd<double, 3> res(s1, nlay, s2, nvar);
  res = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    Range r1(0, pf[j].rows() - 1);
    Range r2(0, pf[j].depth() - 1);
    Array<double, 3> rv(res.value()(r1,ra,r2));
    Array<double, 4> rjac(res.jacobian()(r1,ra,r2,ra));
    Array<double, 1> fv(frac_aer.value()(ra,j));
    Array<double, 2> fjac(frac_aer.jacobian()(ra,j, ra));
    Array<double, 3> pv(pf[j].value());
    Array<double, 4> pjac(pf[j].jacobian());
    rv += fv(i2) * pv(i1, i2, i3);
    if(pf[j].is_constant()) {
      if(!frac_aer.is_constant())
      rjac += fjac(i2,i4) * pv(i1, i2, i3);
    } else {
      if(frac_aer.is_constant())
      rjac += fv(i2) * pjac(i1, i2, i3, i4);
      else
      rjac += fjac(i2,i4) * pv(i1, i2, i3) + fv(i2) * pjac(i1, i2, i3, i4);
    }
  }
  return res;
}

//-----------------------------------------------------------------------
/// Clone a Aerosol object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Pressure>, although you
/// can of course attach them after receiving the cloned object.
//-----------------------------------------------------------------------

boost::shared_ptr<Aerosol> AerosolOptical::clone() const
{
  std::vector<boost::shared_ptr<AerosolExtinction> > aext_clone;
  BOOST_FOREACH(const boost::shared_ptr<AerosolExtinction>& i, aext)
    aext_clone.push_back(i->clone());
  std::vector<boost::shared_ptr<AerosolProperty> > aprop_clone;
  BOOST_FOREACH(const boost::shared_ptr<AerosolProperty>& i, aprop)
    aprop_clone.push_back(i->clone());
  boost::shared_ptr<Aerosol> res(new AerosolOptical(aext_clone, aprop_clone, press->clone(), rh->clone()));
  return res;
}

//-----------------------------------------------------------------------
/// This gives the total aerosol optical depth for a given particle
///
/// You can optionally supply the pressure range to use, this will
/// report the aod for the levels that fall in that range. The default
/// is to use everything.
//-----------------------------------------------------------------------

double AerosolOptical::aerosol_optical_depth(int aer_idx,
                              double pmin, 
                              double pmax) const
{
  double res = 0.0;
  Array<double, 1> pgrid(press->pressure_grid().convert(units::Pa).value.value());
  for(int i = 0; i < pgrid.rows() - 1; ++i) {
    // Three cases. First, everything is in the pmin to pmax range
    if(pmin <= pgrid(i) && pgrid(i + 1) <= pmax)
      res += (pgrid(i + 1) - pgrid(i)) * 
      aext[aer_idx]->extinction_for_layer(i).value();

    // second case is that pmin cuts between pgrid(i) and pgrid(i + 1)
    else if(pmin > pgrid(i) && pmin < pgrid(i + 1) && pgrid(i + 1) <=pmax)
      res += (pgrid(i + 1) - pmin) * aext[aer_idx]->extinction_for_layer(i).value();

    // third case is pmax cuts between pgrid(i) and pgrid(i + 1)
    else if(pmin <= pgrid(i) && pgrid(i) < pmax && pgrid(i + 1) > pmax)
      res += (pmax - pgrid(i)) * aext[aer_idx]->extinction_for_layer(i).value();
  }
  return res;
}

//-----------------------------------------------------------------------
/// This gives the total optical depth for each particle, plus adds
/// the total optical depth for all particles as the 1st entry. This
/// may seem a bit odd, but this is what is expected as output to the
/// HDF file.
///
/// You can optionally supply the pressure range to use, this will
/// report the aod for the levels that fall in that range. The default
/// is to use everything.
//-----------------------------------------------------------------------

double AerosolOptical::aerosol_optical_depth_total
(double pmin, double pmax) const
{
  double res = 0.0;
  for(int j = 0; j < number_particle(); ++j)
    res += aerosol_optical_depth(j, pmin, pmax);
  return res;
}

//-----------------------------------------------------------------------
/// Name of aerosols.
//-----------------------------------------------------------------------
std::vector<std::string> AerosolOptical::aerosol_name() const 
{
  std::vector<std::string> res;
  BOOST_FOREACH(const boost::shared_ptr<AerosolExtinction>& i, aext)
    res.push_back(i->aerosol_name());
  return res;
}

void AerosolOptical::print(std::ostream& Os) const 
{ 
  OstreamPad opad(Os, "    ");
  Os << "AerosolOptical:";
  for(int i = 0; i < (int) aprop.size(); ++i) {
    Os << "  Aerosol Optical Property[" << i << "]:\n";
    opad << *aprop[i] << "\n";
    opad.strict_sync();
    Os << "  Aerosol Extinction[" << i << "]:\n";
    opad << *aext[i] << "\n";
    opad.strict_sync();
  }
}
