#include "aerosol.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Aerosol)
.def(luabind::constructor<const std::vector<boost::shared_ptr<AerosolExtinction> >&,
     const std::vector<boost::shared_ptr<AerosolProperty> >&,
     const boost::shared_ptr<Pressure>&>())
.def("number_particle", &Aerosol::number_particle)
.def("aerosol_extinction", &Aerosol::aerosol_extinction)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for Aerosol
//-----------------------------------------------------------------------

AccumulatedTimer Aerosol::timer("Aerosol");

//-----------------------------------------------------------------------
/// Create an aerosol.
/// \param Aext Aerosol extinction for each aerosol.
/// \param Aerosol_prop Aerosol properties for each aerosol.
/// \param Press The Pressure object that gives the pressure grid.
/// \param Reference_wn The wavenumber that Aext is given for. This
///    is optional, the default value matches the reference band given
///    in the ATB.
//-----------------------------------------------------------------------

Aerosol::Aerosol(const std::vector<boost::shared_ptr<AerosolExtinction> >& Aext,
	 const std::vector<boost::shared_ptr<AerosolProperty> >& Aerosol_prop,
	 const boost::shared_ptr<Pressure>& Press,
		 double Reference_wn)
: aext(Aext),
  aprop(Aerosol_prop),
  press(Press), 
  cache_is_stale(true),
  // Need at least 1 jacobian var such as when not running a retrieval
  nvar(1)
{
  if((int) aprop.size() != number_particle())
    throw Exception("aprop needs to be size of number_particle()");
  qext_at_ref.resize(number_particle());
  for(int i = 0; i < number_particle(); ++i) {
    qext_at_ref(i) = aprop[i]->extinction_coefficient(Reference_wn);
    aext[i]->add_observer(*this);
  }
  press->add_observer(*this);
}

void Aerosol::notify_add(StateVector& Sv)
{
  BOOST_FOREACH(boost::shared_ptr<AerosolExtinction>& i, aext)
    Sv.add_observer(*i);
}

void Aerosol::notify_remove(StateVector& Sv)
{
  BOOST_FOREACH(boost::shared_ptr<AerosolExtinction>& i, aext)
    Sv.remove_observer(*i);
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

blitz::Array<std::string, 1> Aerosol::aerosol_name_arr() const
{
  blitz::Array<std::string, 1> res((int) aext.size());
  for(int i = 0; i < res.rows(); i++)
    res(i) = aext[i]->aerosol_name();
  return res;
}

//-----------------------------------------------------------------------
/// We cache the part of the optical_depth_each_layer calculation that
/// is independent of wn.
//-----------------------------------------------------------------------

void Aerosol::fill_cache() const
{
  if(!cache_is_stale)
    return;
  od_ind_wn.resize(press->number_layer(), number_particle(), nvar);
  for(int i = 0; i < od_ind_wn.rows(); ++i) {
    AutoDerivativeWithUnit<double> delta_press = 
      (press->pressure_grid()(i + 1) - press->pressure_grid()(i));
    AutoDerivative<double> dp = delta_press.convert(units::Pa).value;

    // Resize number of variables in case surface pressure is not retrieved
    if (dp.number_variable() == 0)
      dp.gradient().resize(nvar);

    for(int j = 0; j < number_particle(); ++j) {
      od_ind_wn(i, j) = 1.0 / qext_at_ref(j) * 
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
/// This gives the optical depth for each layer, for the given wave
/// number. Note this only includes the aerosol portion of this,
/// Atmosphere class combines this with Absorbers and rayleigh
/// scattering.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer() x number_particle().
//-----------------------------------------------------------------------

ArrayAd<double, 2> 
Aerosol::optical_depth_each_layer(double wn) const
{
  Range ra(Range::all());
  FunctionTimer ft(timer.function_timer());
  fill_cache();
  ArrayAd<double, 2> res(od_ind_wn.copy());
  for(int i = 0; i < number_particle(); ++i) {
    double t = aprop[i]->extinction_coefficient(wn);
    res.value()(ra, i) *= t;
    res.jacobian()(ra, i, ra) *= t;
  }
  return res;
}

//-----------------------------------------------------------------------
/// This gives the single scatter albedo for each layer, for the given wave
/// number, for the given particle. Note this only includes the
/// aerosol portion of this, 
/// Atmosphere class combines this with Rayleigh scattering.
///
/// We take in the optical depth of each layer. This is just what is
/// returned by optical_depth_each_layer(), we take this in because
/// we can change what the derivative of optical_depth_each_layer is
/// respect to, e.g. in AtmosphereOco we use taua_i.
///
/// This calculates the derivative with respect to whatever variables
/// Od is relative to.
///
/// This has size of number_active_layer()
//-----------------------------------------------------------------------

ArrayAd<double, 1> Aerosol::ssa_each_layer
(double wn,
 int particle_index,
 const ArrayAd<double, 1>& Od) const
{
  FunctionTimer ft(timer.function_timer());
  ArrayAd<double, 1> res(Od.copy());
  double t = aprop[particle_index]->scattering_coefficient(wn) / 
    aprop[particle_index]->extinction_coefficient(wn);
  res.value() *= t;
  if(!res.is_constant())
    res.jacobian() *= t;
  return res;
}

//-----------------------------------------------------------------------
/// This gives the single scatter albedo for each layer, for the given wave
/// number. Note this only includes the
/// aerosol portion of this, 
/// Atmosphere class combines this with Rayleigh scattering.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer()
//-----------------------------------------------------------------------

ArrayAd<double, 1> 
Aerosol::ssa_each_layer(double wn) const
{
  FunctionTimer ft(timer.function_timer());
  ArrayAd<double, 2> od(optical_depth_each_layer(wn));
  ArrayAd<double, 1> res(od.rows(), nvar);
  res.value() = 0;
  res.jacobian() = 0;
  for(int i = 0; i < number_particle(); ++i) {
    ArrayAd<double, 1> t(ssa_each_layer(wn, i, od(Range::all(), i)));
    res.value() += t.value();
    res.jacobian() += t.jacobian();
  }
  return res;
}

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol for a single particle.
/// \param wn The wave number.
/// \param pindex The particle index.
//-----------------------------------------------------------------------

blitz::Array<double, 2> Aerosol::pf_mom(double wn, int pindex) const
{
  FunctionTimer ft(timer.function_timer());
  range_check(pindex, 0, number_particle());
  return aprop[pindex]->phase_function_moment(wn);
}

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol.
/// \param wn The wave number.
/// \param frac_aer This is number_active_layer() x number_particle()
//-----------------------------------------------------------------------

blitz::Array<double, 3> Aerosol::pf_mom(double wn, 
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
  std::vector<Array<double, 2> > pf;
  int s1 = 0;
  int s2 = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    pf.push_back(pf_mom(wn, j));
    s1 = std::max(s1, pf[j].rows());
    s2 = std::max(s2, pf[j].cols());
  }
  blitz::Array<double, 3> res(s1, nlay, s2);
  res = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    Range r1(0, pf[j].rows() - 1);
    Range r2(0, pf[j].cols() - 1);
    res(r1,ra,r2) += frac_aer(ra,j)(i2) * pf[j](i1, i3);
  }
  return res;
}

ArrayAd<double, 3> Aerosol::pf_mom(double wn, 
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

  std::vector<Array<double, 2> > pf;
  int s1 = 0;
  int s2 = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    pf.push_back(aprop[j]->phase_function_moment(wn, nummom, numscat));
    s1 = std::max(s1, pf[j].rows());
    s2 = std::max(s2, pf[j].cols());
  }
  ArrayAd<double, 3> res(s1, nlay, s2, 
			 frac_aer.number_variable());
  res = 0;
  for(int j = 0; j < frac_aer.cols(); ++j) {
    Range r1(0, pf[j].rows() - 1);
    Range r2(0, pf[j].cols() - 1);
    res.value()(r1,ra,r2) += frac_aer.value()(ra,j)(i2) * pf[j](i1, i3);
    res.jacobian()(r1,ra,r2,ra) += 
      frac_aer.jacobian()(ra,j,ra)(i2,i4) * pf[j](i1, i3);
  }
  return res;
}

//-----------------------------------------------------------------------
/// This version of clone takes a pressure to use. The intent is that
/// the pressure has been cloned from the original pressure (although
/// this class has no way to verify this). This allows sets of objects
/// to be cloned using a common Pressure clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

boost::shared_ptr<Aerosol> 
Aerosol::clone(const boost::shared_ptr<Pressure>& Press) const
{
  std::vector<boost::shared_ptr<AerosolExtinction> > aext_clone;
  BOOST_FOREACH(const boost::shared_ptr<AerosolExtinction>& i, aext)
    aext_clone.push_back(i->clone(Press));
  boost::shared_ptr<Aerosol> res(new Aerosol(aext_clone, aprop, Press));
  return res;
}

//-----------------------------------------------------------------------
/// This gives the total aerosol optical depth for a given particle
///
/// You can optionally supply the pressure range to use, this will
/// report the aod for the levels that fall in that range. The default
/// is to use everything.
//-----------------------------------------------------------------------

double Aerosol::aerosol_optical_depth(int aer_idx,
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

double Aerosol::aerosol_optical_depth_total
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
std::vector<std::string> Aerosol::aerosol_name() const 
{
  std::vector<std::string> res;
  BOOST_FOREACH(const boost::shared_ptr<AerosolExtinction>& i, aext)
    res.push_back(i->aerosol_name());
  return res;
}

void Aerosol::print(std::ostream& Os) const 
{ 
  OstreamPad opad(Os, "    ");
  Os << "Aerosol:";
  for(int i = 0; i < (int) aprop.size(); ++i) {
    Os << "  Aerosol Property[" << i << "]:\n";
    opad << *aprop[i] << "\n";
    opad.strict_sync();
    Os << "  Aerosol Extinction[" << i << "]:\n";
    opad << *aext[i] << "\n";
    opad.strict_sync();
  }
}
