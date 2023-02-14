#include "vlidort_rt_driver.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ground.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void VLidortRtDriver::serialize(Archive & ar, const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiScattRtDriver)
    & FP_NVP_(vlidort_interface);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidortRtDriver::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void VLidortRtDriver::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // Recreate RT pars, no need to serialize this object
  rt_pars_.reset( new VLidort_Pars() );
}

FP_IMPLEMENT(VLidortRtDriver);

#endif

//=======================================================================
// VLidortRtDriver
//=======================================================================

VLidortRtDriver::VLidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, 
                                 int surface_type, const blitz::Array<double, 1>& zen, bool pure_nadir,
                                 bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering)
  : MultiScattRtDriver(do_solar_sources, do_thermal_emission)
{
  initialize_interface(nstream, nmoment);
  initialize_brdf(surface_type);
  initialize_rt(nstream, nmoment, do_solar_sources, do_thermal_emission, do_thermal_scattering);
  
  // Set up scatting mode based on viewing zenith angle
  setup_sphericity(zen, do_multi_scatt_only, pure_nadir);
}

void VLidortRtDriver::initialize_interface(int nstream, int nmoment)
{
  rt_pars_.reset( new VLidort_Pars() );
  brdf_driver_.reset( new VLidortBrdfDriver(nstream, nmoment) );
  vlidort_interface_.reset( new VLidort_Lps_Masters() );

  // Check inputs against sizes allowed by LIDORT
  range_check(nstream, 1, rt_pars_->maxstreams()+1);
  range_check(nmoment, 2, rt_pars_->maxmoments_input()+1);
}

void VLidortRtDriver::setup_phase_function(const blitz::Array<double, 3>& pf)
{
  // Ranges for copying inputs to method
  Range rall = Range::all();
  Range rlay(0, pf.extent(secondDim) - 1);
  Range rmom(0, pf.extent(firstDim) - 1);

  /*
  // For all layers n and threads t, Legrenre moments of
  // the phase function expansion multiplied by (2L+1);
  // initial value (L=0) should always be 1
  // phasmoms_total_input(n, L, t)
  // n = moments, L = layers, t = threads
  Lidort_Fixed_Optical& lid_foptical_inputs = vlidort_interface_->lidort_fixin().optical();

  Array<double, 2> pf_scalar(pf(rall, rall, 0));

  Array<double, 2> phasmoms( lid_foptical_inputs.ts_phasmoms_total_input() );
  phasmoms(rmom, rlay) = where(abs(pf_scalar) > 1e-11, pf_scalar, 1e-11);
  */
}

void VLidortRtDriver::setup_linear_phase_function(const ArrayAd<double, 3>& pf)
{

  int natm_jac = pf.number_variable();

  // Ranges for copying inputs to method
  Range rall = Range::all();
  Range rmom(0, pf.rows() - 1); // number phase function moments
  Range rlay(0, pf.cols() - 1);
  Range rjac(0, natm_jac - 1);
  firstIndex i1; secondIndex i2;

/*
  Lidort_Fixed_Linoptical& lid_linoptical = vlidort_interface_->lidort_linfixin().optical();
  Array<double, 3> l_phasmoms( lid_linoptical.ts_l_phasmoms_total_input()(rjac,rmom,rlay) );

  if(pf.is_constant())
    l_phasmoms(rjac, rmom, rlay) = 0.0;
  else {
    Array<double, 2> pf_val_scalar(pf.value()(rall, rall, 0));
    blitz::Array<double, 2> pf_in( where(abs(pf_val_scalar) > 1e-11, pf_val_scalar, 1e-11) );

    // We need this loop since l_phasmoms and pf variables have jacobian data in different dimensions
    Array<double, 3> pf_jac_scalar(pf.jacobian()(rall, rall, 0, rall));
    for (int jidx = 0; jidx < pf.number_variable(); jidx++)
      l_phasmoms(jidx, rmom, rlay) = pf_jac_scalar(rmom, rlay, jidx) / pf_in(rmom, rlay)(i1, i2);
  }
*/
}

/// Copy outputs from BRDF supplement into LIDORT Sup inputs types
void VLidortRtDriver::copy_brdf_sup_outputs() const {
/*
  // Copy BRDF outputs to LIDORT's BRDF inputs
  // Use LIDORT specific interfaces to ensure copying is done in memory correctly
  // Cannot mix VLIDORT and LIDORT BRDF objects here
  Brdf_Sup_Outputs& brdf_outputs = lidort_brdf_interface()->brdf_sup_out();
  Brdf_Linsup_Outputs& brdf_lin_outputs = lidort_brdf_interface()->brdf_linsup_out();

  Lidort_Sup_Brdf& lid_brdf = vlidort_interface_->lidort_sup().brdf();
  Lidort_Linsup_Brdf& lid_lin_brdf = vlidort_interface_->lidort_linsup().brdf();

  lid_brdf.copy_from_sup( brdf_outputs );
  lid_lin_brdf.copy_from_sup( brdf_lin_outputs );
*/
}

const blitz::Array<double, 1> VLidortRtDriver::get_intensity() const
{
/*
  // Total Intensity I(t,v,d,T) at output level t, output geometry v,
  // direction d
  Array<double, 1> intensity(1);
  intensity(0) = vlidort_interface_->lidort_out().main().ts_intensity()(0,0, rt_pars_->upidx()-1);
  return intensity;
*/
}

void VLidortRtDriver::copy_jacobians(blitz::Array<double, 3>& jac_atm, blitz::Array<double, 2>& jac_surf_param, blitz::Array<double, 1>& jac_surf_temp, blitz::Array<double, 2>& jac_atm_temp) const
{
/*
  Lidort_Linatmos& lpoutputs = vlidort_interface_->lidort_linout().atmos();
  Lidort_Linsurf& lsoutputs = vlidort_interface_->lidort_linout().surf();

  Range ra(Range::all());

  // Surface Jacobians KR(r,t,v,d) with respect to surface variable r
  // at output level t, geometry v, direction d
  Array<double, 1> jac_surf_param_interface( lsoutputs.ts_surfacewf()(ra, 0, 0, rt_pars_->upidx()-1) );
  jac_surf_param.resize(jac_surf_param_interface.rows(), 1);
  jac_surf_param(ra, 0) = jac_surf_param_interface;

  // Get profile jacobians
  // Jacobians K(q,n,t,v,d) with respect to profile atmospheric variable
  // q in layer n, at output level t, geometry v, direction d
  Array<double, 2> jac_atm_interface( lpoutputs.ts_profilewf()(ra, ra, 0, 0, rt_pars_->upidx()-1) );
  jac_atm.resize(jac_atm_interface.rows(), jac_atm_interface.cols(), 1);
  jac_atm(ra, ra, 0) = jac_atm_interface;

  // Get surface temp jacobian if thermal emission is enabled
  if(do_thermal_emission) {
      jac_surf_temp.resize(1);
      jac_surf_temp(0) = lsoutputs.ts_sbbwfs_jacobians()(0, 0, rt_pars_->upidx()-1);

      Array<double, 1> jac_atm_temp_interface( lpoutputs.ts_abbwfs_jacobians()(0, 0, ra, rt_pars_->upidx()-1) );
      jac_atm_temp.resize(jac_atm_temp_interface.rows(), 1);
      jac_atm_temp(ra, 0) = jac_atm_temp_interface;
  }
*/
}
