#include "chapman_boa_rt.h"
#include "fp_serialize_support.h"
#include "wgs84_constant.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ChapmanBoaRTCache::serialize(Archive & ar,
				  const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CacheInvalidatedObserver)
    & FP_NVP(chapman_boa);
}

template<class Archive>
void ChapmanBoaRT::serialize(Archive & ar,
			     const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RadiativeTransfer)
    &  FP_NVP(spec_bound) & FP_NVP(atm) & FP_NVP(sza);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void ChapmanBoaRT::save(Archive &UNUSED(ar),
		       const unsigned int UNUSED(version)) const
{
}

template<class Archive>
void ChapmanBoaRT::load(Archive &UNUSED(ar),
			const unsigned int UNUSED(version)) 
{
  init();
}

FP_IMPLEMENT(ChapmanBoaRTCache);
FP_IMPLEMENT(ChapmanBoaRT);
#endif

boost::shared_ptr<RadiativeTransfer> chapman_boa_rt_create
(const boost::shared_ptr<RtAtmosphere>& Atm,
 const blitz::Array<double, 1>& Sza, 
 const SpectralBound& Spec_bound)
{
  const boost::shared_ptr<AtmosphereStandard> atm_oco(boost::dynamic_pointer_cast<AtmosphereStandard>(Atm));
  return boost::shared_ptr<RadiativeTransfer>(new ChapmanBoaRT(atm_oco, Sza, Spec_bound));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ChapmanBoaRT, RadiativeTransfer)
.scope
[
 luabind::def("create", &chapman_boa_rt_create)
]
REGISTER_LUA_END()
#endif

ChapmanBoaRT::ChapmanBoaRT(const boost::shared_ptr<AtmosphereStandard>& Atm,
			   const blitz::Array<double, 1>& Sza) 
  : atm(Atm), sza(Sza)
{
  init();
}

void ChapmanBoaRT::init()
{
  cache.resize(sza.rows());
  BOOST_FOREACH(auto i, cache)
    atm->add_cache_invalidated_observer(i);
}

ChapmanBoaRT::ChapmanBoaRT(const boost::shared_ptr<AtmosphereStandard>& Atm,
			   const blitz::Array<double, 1>& Sza, 
			   const SpectralBound& Spec_bound) 
  : spec_bound(Spec_bound), atm(Atm), sza(Sza)
{
  cache.resize(spec_bound.number_spectrometer());
  BOOST_FOREACH(auto i, cache)
    atm->add_cache_invalidated_observer(i);
}

void ChapmanBoaRTCache::fill_cache(const ChapmanBoaRT& C, int Spec_index)
{
  double rearth = OldConstant::wgs84_a.convert(units::km).value;
  double rfindex_param = 0.000288;

  // Can not use OCO refracrive index if CO2 and H2O are not
  // present in the state structure
  bool can_use_oco_refr = C.atm->absorber_ptr()->gas_index("CO2") != -1 and 
    C.atm->absorber_ptr()->gas_index("H2O") != -1;

  // Atmospheric values
  Array<AutoDerivative<double>, 1> height_grid( C.atm->altitude(Spec_index).convert(units::km).value.to_array() );
  Array<AutoDerivative<double>, 1> press_grid( C.atm->pressure_ptr()->pressure_grid().convert(units::Pa).value.to_array() );
  Array<AutoDerivative<double>, 1> temp_grid(press_grid.rows());
  Array<AutoDerivative<double>, 1> co2_vmr(press_grid.rows());
  Array<AutoDerivative<double>, 1> h2o_vmr(press_grid.rows());
  for(int i = 0; i < temp_grid.rows(); ++i) {
    temp_grid(i) =
      C.atm->temperature_ptr()->temperature(AutoDerivativeWithUnit<double>(press_grid(i), units::Pa)).convert(units::K).value;
    
    if(can_use_oco_refr) {
      co2_vmr(i) = C.atm->absorber_ptr()->absorber_vmr("CO2")->
	volume_mixing_ratio(press_grid(i));
      h2o_vmr(i) = C.atm->absorber_ptr()->absorber_vmr("H2O")->
	volume_mixing_ratio(press_grid(i));
    }
  }
  
  // Calculate reference wavelengths, if we can..
  double ref_wavelength = -1;
  if(C.spec_bound.number_spectrometer() > 0)
    ref_wavelength = C.spec_bound.center(Spec_index, units::micron).value;
  
  Logger::debug() << "Generating Chapman Factors\n";
  boost::shared_ptr<AtmRefractiveIndex> refr_index;
  Logger::debug() << "Band " << (Spec_index + 1) << " using ";
  if(can_use_oco_refr && C.spec_bound.number_spectrometer() > 0 && 
     OcoRefractiveIndex::wavelength_in_bounds(ref_wavelength)) {
    Logger::debug() << "OCO";
    refr_index.reset(new OcoRefractiveIndex(ref_wavelength, press_grid,
					    temp_grid, co2_vmr, h2o_vmr));
  } else {
    Logger::debug() << "simple";
    refr_index.reset(new SimpleRefractiveIndex(rfindex_param,
					       press_grid, temp_grid));
  }
  Logger::debug() << " refractive index class.\n";
  
  chapman_boa = boost::make_shared<ChapmanBOA>
    (rearth, C.sza(Spec_index), height_grid, refr_index); 
}

Spectrum ChapmanBoaRT::reflectance(const SpectralDomain& Spec_domain, 
				int Spec_index, bool Skip_jacobian) const
{
  ArrayAd<double, 1> res;
  if(Skip_jacobian)
    res.reference(ArrayAd<double, 1>(stokes(Spec_domain, Spec_index)
				     (Range::all(), 0)));
  else 
    res.reference(stokes_and_jacobian(Spec_domain, Spec_index)
		  (Range::all(), 0));
  return Spectrum(Spec_domain, SpectralRange(res, units::inv_sr));
}

blitz::Array<double, 2> ChapmanBoaRT::stokes(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));
  /// \todo : Make this faster to where it doesn't compute jacobians, right now it just
  /// fails to return the jacobian part and it is computed unncecessarily

  // Compute chapman factors if needed
  cache[Spec_index].fill_cache_if_needed(*this, Spec_index);

  Array<double, 1> wn(Spec_domain.wavenumber());
  Array<double, 2> trans(wn.extent(firstDim), number_stokes());
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);

  for(int i = 0; i < wn.rows(); ++i) {
    trans(i, 0) = value( cache[Spec_index].chapman_boa->transmittance(atm->optical_depth_wrt_state_vector(wn(i), Spec_index).to_array(), 0) );
    if(disp) *disp += 1;
  }

  return trans;
}

ArrayAd<double, 2> ChapmanBoaRT::stokes_and_jacobian(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));

  // Compute chapman factors if needed
  cache[Spec_index].fill_cache_if_needed(*this, Spec_index);

  Array<double, 1> wn(Spec_domain.wavenumber());
  ArrayAd<double, 2> trans_jac(wn.extent(firstDim), number_stokes(), 1);
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);

  for(int i = 0; i < wn.extent(firstDim); ++i) {
    AutoDerivative<double> trans_wn = cache[Spec_index].chapman_boa->transmittance(atm->optical_depth_wrt_state_vector(wn(i), Spec_index).to_array(), 0);
    trans_jac.resize_number_variable(trans_wn.number_variable());
    trans_jac(i, 0) = trans_wn;
    if(disp) *disp += 1;
  }

  return trans_jac;
}

void ChapmanBoaRT::print(std::ostream& Os, bool Short_form) const
{
  Os << "ChapmanBoaRT";
  OstreamPad opad(Os, "  ");
  if(!Short_form) {
    Os << "\nAtmosphere:\n";
    opad << *atm << "\n";
    opad.strict_sync();
  }
}
