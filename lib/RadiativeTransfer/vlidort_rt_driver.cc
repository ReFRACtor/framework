#include "vlidort_rt_driver.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ground.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
// Mapping of the 6 de-Rooij scattering-matrix expansion coefficients
// carried in the third dimension of the phase function array (order
// [a1, a2, a3, a4, b1, b2], see RtAtmosphere::phase_function_moments_wrt_rt
// and RayleighGreekMoment) into the 16-element (MAXSTOKES_SQ) VLIDORT
// "greek matrix" input (ts_greekmat_total_input). This reproduces the
// GMASK/RMASK masking used by the VLIDORT 2.8.3 tester
// (test_vlidort_lin_2p8p3.f90):
//
//   GMASK = { 1, 2, 5, 6, 11, 12, 15, 16 }   ! 1-based greek slot
//   RMASK = { 1, 5, 5, 2,  3,  6,  6,  4 }   ! 1-based de-Rooij element
//
// The number of entries used depends on nstokes: 1 (nstokes=1),
// 5 (nstokes=3) and 8 (nstokes=4).
//
// Note on sign: RefractOr's b1 (gamma coupling) coefficient uses the
// opposite sign convention from the value VLIDORT expects in the (1,2)/
// (2,1) greek slots (RayleighGreekMoment uses +sqrt(6)*beta2 while VLIDORT
// wants -sqrt(6)*beta2). The b1 element is therefore negated when copied
// into the greek slots so that the Stokes output matches the VLIDORT
// benchmark tester.
//-----------------------------------------------------------------------

namespace {
  // 0-based greek slot for each masked entry
  const int greek_slot[8] = { 0, 1, 4, 5, 10, 11, 14, 15 };
  // 0-based de-Rooij element feeding each masked entry
  const int derooij_elem[8] = { 0, 4, 4, 1, 2, 5, 5, 3 };

  int num_greek_entries(int nstokes)
  {
    if (nstokes == 1) return 1;
    if (nstokes == 3) return 5;
    if (nstokes == 4) return 8;
    throw FullPhysics::Exception("VLIDORT nstokes must be 1, 3 or 4");
  }
}

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void VLidortRtDriver::serialize(Archive & ar, const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiScattRtDriver)
    & FP_NVP_(vlidort_interface) & FP_NVP_(nstokes);
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
                                 int nstokes,
                                 bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering)
  : MultiScattRtDriver(do_solar_sources, do_thermal_emission), nstokes_(nstokes)
{
  range_check(nstokes, 1, 5);
  if(nstokes == 2)
    throw Exception("VLIDORT nstokes must be 1, 3 or 4");

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

  // Check inputs against sizes allowed by VLIDORT
  range_check(nstream, 1, rt_pars_->maxstreams()+1);
  range_check(nmoment, 2, rt_pars_->maxmoments_input()+1);
}

void VLidortRtDriver::initialize_rt(int nstream, int nmoment, bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering)
{
  MultiScattRtDriver::initialize_rt(nstream, nmoment, do_solar_sources, do_thermal_emission, do_thermal_scattering);

  // VLIDORT specific
  VLidort_Fixed_Control& fcontrol_inputs = vlidort_interface_->vlidort_fixin().cont();
  fcontrol_inputs.ts_nstokes(nstokes_);

  // The VBRDF supplement needs to know the number of Stokes parameters so
  // it computes the full polarized surface reflectance matrix. Without this
  // it defaults to zero and no surface contribution is produced.
  vlidort_brdf_interface()->vbrdf_sup_in().bs_nstokes(nstokes_);

  // Always use the BRDF supplement, don't use the specialized
  // lambertian_albedo mode (matches the LIDORT driver's use of
  // ts_do_brdf_surface(true))
  VLidort_Fixed_Boolean& fboolean_inputs = vlidort_interface_->vlidort_fixin().f_bool();
  fboolean_inputs.ts_do_lambertian_surface(false);

  // Special Fourier convergence acceleration when nstokes == 2, harmless
  // otherwise. Matches the VLIDORT 2.8.3 tester configuration.
  fboolean_inputs.ts_do_fourier0_nstokes2(true);
}

void VLidortRtDriver::setup_phase_function(const blitz::Array<double, 3>& pf)
{
  // Ranges for copying inputs to method
  Range rlay(0, pf.extent(secondDim) - 1);
  Range rmom(0, pf.extent(firstDim) - 1);

  // For all layers n and moments L, expansion coefficients of the greek
  // scattering matrix. ts_greekmat_total_input(L, n, k), k = 0..15 indexes
  // the flattened 4x4 greek matrix. The (1,1) element (slot 0) is the
  // traditional phase function moment used by the scalar RTE.
  VLidort_Fixed_Optical& foptical_inputs = vlidort_interface_->vlidort_fixin().optical();

  Array<double, 3> greekmat( foptical_inputs.ts_phasmoms_total_input() );

  int n_entries = num_greek_entries(nstokes_);
  for(int k = 0; k < n_entries; ++k) {
    int gk = greek_slot[k];
    int rk = derooij_elem[k];
    Array<double, 2> src( pf(rmom, rlay, rk) );
    // b1 (de-Rooij element index 4) uses the opposite sign convention from
    // the VLIDORT greek off-diagonal (1,2)/(2,1) slots.
    double sign = (rk == 4) ? -1.0 : 1.0;
    if(gk == 0)
      // (1,1) element: floor small/zero values so VLIDORT does not choke,
      // exactly as the LIDORT driver does for phasmoms_total_input.
      greekmat(rmom, rlay, gk) = where(abs(src) > 1e-11, src, 1e-11);
    else
      greekmat(rmom, rlay, gk) = sign * src;
  }
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

  VLidort_Fixed_Linoptical& lin_optical = vlidort_interface_->vlidort_linfixin().optical();
  Array<double, 4> l_greekmat( lin_optical.ts_l_phasmoms_total_input() );

  if(pf.is_constant()) {
    int n_entries = num_greek_entries(nstokes_);
    for(int k = 0; k < n_entries; ++k)
      l_greekmat(rjac, rmom, rlay, greek_slot[k]) = 0.0;
    return;
  }

  int n_entries = num_greek_entries(nstokes_);
  for(int k = 0; k < n_entries; ++k) {
    int gk = greek_slot[k];
    int rk = derooij_elem[k];
    double sign = (rk == 4) ? -1.0 : 1.0;

    Array<double, 2> pf_val( pf.value()(rall, rall, rk) );
    // Floor the (1,1) element consistently with setup_phase_function so the
    // relative-derivative division below is well defined.
    Array<double, 2> pf_in( (gk == 0) ? Array<double,2>(where(abs(pf_val) > 1e-11, pf_val, 1e-11))
                                      : pf_val );

    // l_greekmat and pf carry jacobian data in different dimensions, so loop
    // over the jacobian index. VLIDORT stores the linearization as the
    // relative derivative (dGreek/dx)/Greek.
    Array<double, 3> pf_jac( pf.jacobian()(rall, rall, rk, rall) );
    for (int jidx = 0; jidx < natm_jac; ++jidx) {
      if(gk == 0)
        l_greekmat(jidx, rmom, rlay, gk) = pf_jac(rmom, rlay, jidx) / pf_in(rmom, rlay)(i1, i2);
      else
        l_greekmat(jidx, rmom, rlay, gk) = sign * pf_jac(rmom, rlay, jidx) / pf_in(rmom, rlay)(i1, i2);
    }
  }
}

/// Copy outputs from BRDF supplement into VLIDORT Sup inputs types
void VLidortRtDriver::copy_brdf_sup_outputs() const {

  // Copy BRDF outputs to VLIDORT's BRDF inputs
  // Use VLIDORT specific interfaces to ensure copying is done in memory correctly
  // Cannot mix VLIDORT and LIDORT BRDF objects here
  VBrdf_Sup_Outputs& brdf_outputs = vlidort_brdf_interface()->vbrdf_sup_out();
  VBrdf_Linsup_Outputs& brdf_lin_outputs = vlidort_brdf_interface()->vbrdf_linsup_out();

  VLidort_Sup_Brdf& vlid_brdf = vlidort_interface_->vlidort_sup().brdf();
  VLidort_Linsup_Brdf& vlid_lin_brdf = vlidort_interface_->vlidort_linsup().brdf();

  vlid_brdf.copy_from_sup( brdf_outputs );
  vlid_lin_brdf.copy_from_sup( brdf_lin_outputs );
}

const blitz::Array<double, 1> VLidortRtDriver::get_intensity() const
{
  // Stokes vector I (and Q, U, V) at output level t, output geometry v,
  // upwelling direction. ts_stokes(t, v, stokes, direction).
  Array<double, 1> stokes(nstokes_);
  const Array<double, 4>& ts_stokes = vlidort_interface_->vlidort_out().main().ts_stokes();
  for(int s = 0; s < nstokes_; ++s)
    stokes(s) = ts_stokes(0, 0, s, rt_pars_->upidx()-1);
  return stokes;
}

void VLidortRtDriver::copy_jacobians(blitz::Array<double, 3>& jac_atm, blitz::Array<double, 2>& jac_surf_param, blitz::Array<double, 1>& jac_surf_temp, blitz::Array<double, 2>& jac_atm_temp) const
{
  VLidort_Linprof& lpoutputs = vlidort_interface_->vlidort_linout().prof();
  VLidort_Linatmos& latmos = vlidort_interface_->vlidort_linout().atmos();
  VLidort_Linsurf& lsoutputs = vlidort_interface_->vlidort_linout().surf();

  Range ra(Range::all());
  int upidx = rt_pars_->upidx()-1;

  // Surface Jacobians KR(r,t,v,stokes,d) with respect to surface variable r
  // at output level t, geometry v, direction d
  const Array<double, 5>& surfacewf = lsoutputs.ts_surfacewf();
  int n_surf = surfacewf.extent(firstDim);
  jac_surf_param.resize(n_surf, nstokes_);
  for(int s = 0; s < nstokes_; ++s)
    jac_surf_param(ra, s) = surfacewf(ra, 0, 0, s, upidx);

  // Profile jacobians K(q,n,t,v,stokes,d) with respect to profile
  // atmospheric variable q in layer n, at output level t, geometry v,
  // direction d
  const Array<double, 6>& profilewf = lpoutputs.ts_profilewf();
  int n_jac = profilewf.extent(firstDim);
  int n_lay = profilewf.extent(secondDim);
  jac_atm.resize(n_jac, n_lay, nstokes_);
  for(int s = 0; s < nstokes_; ++s)
    jac_atm(ra, ra, s) = profilewf(ra, ra, 0, 0, s, upidx);

  // Get thermal temperature jacobians if thermal emission is enabled
  if(do_thermal_emission) {
      const Array<double, 4>& sbbwfs = lsoutputs.ts_sbbwfs_jacobians();
      const Array<double, 5>& abbwfs = latmos.ts_abbwfs_jacobians();

      jac_surf_temp.resize(nstokes_);
      jac_atm_temp.resize(abbwfs.extent(thirdDim), nstokes_);
      for(int s = 0; s < nstokes_; ++s) {
          jac_surf_temp(s) = sbbwfs(0, 0, s, upidx);
          jac_atm_temp(ra, s) = abbwfs(0, 0, ra, s, upidx);
      }
  }
}
