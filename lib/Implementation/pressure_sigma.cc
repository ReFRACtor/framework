#include "pressure_sigma.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void PressureSigma::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PressureImpBase)
    & FP_NVP_(a) & FP_NVP_(b);
}

FP_IMPLEMENT(PressureSigma);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(PressureSigma, Pressure)
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<double, 1>&,
                          double>())
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          double>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

PressureSigma::PressureSigma(const blitz::Array<double, 1>& A,
                             const blitz::Array<double, 1>& B,
                             double Surface_pressure,
			     Pressure::TypePreference Tpref)
: a_(A.copy()), b_(B.copy())
{
  if(A.rows() != B.rows())
    throw Exception("A and B need to be the same size in PressureSigma constructor");

  type_preference_ = Tpref;
  blitz::Array<double, 1> val(1);
  val(0) = Surface_pressure;
  init(val);
  cov.resize(1, 1);
  cov(0,0) = 1;
}

//-----------------------------------------------------------------------
/// Constructor.
/// Creates A and B parameters from the pressure grid passed in. 
/// A becomes all 0 of the same size as Pressure_grid
/// B becomes Pressure_grid / Pressure_grid[-1]
//-----------------------------------------------------------------------

PressureSigma::PressureSigma(const blitz::Array<double, 1>& Pressure_grid,
                             double Surface_pressure,
			     Pressure::TypePreference Tpref)
{
  type_preference_ = Tpref;
  set_levels_from_grid(Pressure_grid);

  blitz::Array<double, 1> val(1);
  val(0) = Surface_pressure;
  init(val);
  cov.resize(1, 1);
  cov(0,0) = 1;
}

//-----------------------------------------------------------------------
/// Creates A and B parameters from the pressure grid passed in. 
/// A becomes all 0 of the same size as Pressure_grid
/// B becomes Pressure_grid / Pressure_grid[-1]
//-----------------------------------------------------------------------

void PressureSigma::set_levels_from_grid(const blitz::Array<double, 1>& Pressure_grid)
{
  a_.resize(Pressure_grid.rows());
  a_ = 0.0;

  b_.resize(Pressure_grid.rows());
  b_ = Pressure_grid;
  if(type_preference_ == Pressure::PREFER_INCREASING_PRESSURE)
    b_ = b_ / Pressure_grid(Pressure_grid.rows()-1);
  else
    b_ = b_ / Pressure_grid(0);
  
  Observable<Pressure>::notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// Calculate the new pressure grid. Done each time the surface
/// pressure is updated.
//-----------------------------------------------------------------------

void PressureSigma::calc_pressure_grid() const
{
  cache.pgrid.units = units::Pa;
  cache.pgrid.value.resize(b_.rows(), coeff.number_variable());
  for(int i = 0; i < b_.rows(); ++i) {
    cache.pgrid.value(i) = b_(i) * coeff(0) + a_(i);

    // Ensure that the pressure grid calculated in increasing order if
    // that is what we said
    // Since so many linear interpolations rely on this, this error
    // message will make more sense than what would be throw otherwise:
    // X needs to be sorted
    if(type_preference() == Pressure::PREFER_INCREASING_PRESSURE &&
       i > 0 &&
       cache.pgrid.value(i-1).value() > cache.pgrid.value(i).value()) {
      stringstream err_msg;
      err_msg << "At level " << i << " pressure is smaller: " 
	      << "(" 
	      << b_(i) << " * " << coeff(0).value() << " + " << a_(i)
	      << ") = "
	      << cache.pgrid.value(i).value()
	      << " than the value at the previous level: " 
	      << "(" 
	      << b_(i-1) << " * " << coeff(0).value() << " + " << a_(i-1)
	      << ") = "
	      << cache.pgrid.value(i-1).value();
      throw Exception(err_msg.str());
    }
    if(type_preference() == Pressure::PREFER_DECREASING_PRESSURE &&
       i > 0 &&
       cache.pgrid.value(i-1).value() < cache.pgrid.value(i).value()) {
      stringstream err_msg;
      err_msg << "At level " << i << " pressure is larger: " 
	      << "(" 
	      << b_(i) << " * " << coeff(0).value() << " + " << a_(i)
	      << ") = "
	      << cache.pgrid.value(i).value()
	      << " than the value at the previous level: " 
	      << "(" 
	      << b_(i-1) << " * " << coeff(0).value() << " + " << a_(i-1)
	      << ") = "
	      << cache.pgrid.value(i-1).value();
      throw Exception(err_msg.str());
    }
  }
}

//-----------------------------------------------------------------------
/// Clone a PressureSigma object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<PressureSigma>, although you
/// can of course attach them after receiving the cloned object.
//-----------------------------------------------------------------------

boost::shared_ptr<Pressure> PressureSigma::clone() const
{
  boost::shared_ptr<Pressure> res
    (new PressureSigma(a_, b_, coefficient()(0).value(), type_preference_));
  return res;
}

void PressureSigma::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "PressureSigma:\n"
     << "  Surface pressure:    " << surface_pressure().value.value() << "\n"
     << "  Type preference:     "
     << (type_preference() == Pressure::PREFER_INCREASING_PRESSURE ?
	 "Increasing pressure" : "Decreasing pressure") << "\n"
     << "  a: \n";
  opad << a() << "\n";
  opad.strict_sync();
  Os << "  b: \n";
  opad << b() << "\n";
  opad.strict_sync();
}

