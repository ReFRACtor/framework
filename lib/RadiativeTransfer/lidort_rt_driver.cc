#include "lidort_rt_driver.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ground.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void LidortRtDriver::serialize(Archive & ar, const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiScattRtDriver)
    & FP_NVP_(do_thermal_scattering)
    & FP_NVP_(lidort_interface);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void LidortRtDriver::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void LidortRtDriver::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // Recreate RT pars, no need to serialize this object
  rt_pars_.reset( new Lidort_Pars() );
}

FP_IMPLEMENT(LidortRtDriver);

#endif

//=======================================================================
// LidortRtDriver
//=======================================================================

LidortRtDriver::LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, 
        int surface_type, const blitz::Array<double, 1>& zen, bool pure_nadir,
        bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering)
  : MultiScattRtDriver(do_solar_sources, do_thermal_emission),
    do_thermal_scattering_(do_thermal_scattering)
{
  initialize_interface(nstream, nmoment);
  initialize_brdf(surface_type);
  initialize_rt(nstream, nmoment, do_solar_sources, do_thermal_emission, do_thermal_scattering);
  
  // Set up scatting mode based on viewing zenith angle
  setup_sphericity(zen, do_multi_scatt_only, pure_nadir);
 
}

void LidortRtDriver::initialize_interface(int nstream, int nmoment)
{
  rt_pars_.reset( new Lidort_Pars() );
  brdf_driver_.reset( new LidortBrdfDriver(nstream, nmoment) );
  lidort_interface_.reset( new Lidort_Lps_Masters() );

  // Check inputs against sizes allowed by LIDORT
  range_check(nstream, 1, rt_pars_->maxstreams()+1);
  range_check(nmoment, 2, rt_pars_->maxmoments_input()+1);

}

void LidortRtDriver::initialize_rt(int nstream, int nmoment, bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering)
{
  MultiScattRtDriver::initialize_rt(nstream, nmoment, do_solar_sources, do_thermal_emission, do_thermal_scattering);
  
  // LIDORT specific
  // Always use BRDF supplement, don't use specialized lambertian_albedo mode
  Lidort_Fixed_Boolean& lid_fboolean_inputs = lidort_interface_->lidort_fixin().f_bool();
  lid_fboolean_inputs.ts_do_brdf_surface(true);

}

/// Set plane parallel sphericity
void LidortRtDriver::set_plane_parallel()
{
  MultiScattRtDriver::set_plane_parallel();

  // LIDORT specific
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  mboolean_inputs.ts_do_no_azimuth(true);
}

/// Set pseudo spherical sphericity
void LidortRtDriver::set_pseudo_spherical()
{
  MultiScattRtDriver::set_pseudo_spherical();

  // LIDORT specific
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  mboolean_inputs.ts_do_no_azimuth(false);
}

/// Set plane parallel plus single scattering correction
void LidortRtDriver::set_plane_parallel_plus_ss_correction()
{
  MultiScattRtDriver::set_plane_parallel_plus_ss_correction();

  // LIDORT specific
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  mboolean_inputs.ts_do_no_azimuth(false);
}

/// Set line of sight mode
void LidortRtDriver::set_line_of_sight()
{
  MultiScattRtDriver::set_line_of_sight();

  // LIDORT specific
  Lidort_Modified_Boolean& mboolean_inputs = lidort_interface_->lidort_modin().mbool();
  mboolean_inputs.ts_do_no_azimuth(false);
}

void LidortRtDriver::setup_height_grid(const blitz::Array<double, 1>& in_height_grid)
{
  Lidort_Fixed_Chapman& fchapman_inputs = lidort_interface_->lidort_fixin().chapman();
  Lidort_Modified_Uservalues& muser_inputs = lidort_interface_->lidort_modin().muserval();

  Array<double, 1> lidort_height_grid( fchapman_inputs.ts_height_grid() );
  int nlayer = in_height_grid.extent(firstDim) - 1;
  Range lay_range = Range(0,nlayer);
  lidort_height_grid(lay_range) = in_height_grid;

  // Set GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
  muser_inputs.ts_geometry_specheight( lidort_height_grid(nlayer) );

  // Tell LIDORT number of layers
  lidort_interface_->lidort_fixin().cont().ts_nlayers(nlayer);
}

void LidortRtDriver::setup_geometry(double sza, double azm, double zen)
{

  Lidort_Modified_Sunrays& mbeam_inputs = lidort_interface_->lidort_modin().msunrays();
  Lidort_Modified_Uservalues& muser_inputs = lidort_interface_->lidort_modin().muserval();

  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> ld_sza( mbeam_inputs.ts_beam_szas() );
  ld_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> ld_azm( muser_inputs.ts_user_relazms() );
  ld_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> ld_zen( muser_inputs.ts_user_angles_input() );
  ld_zen(0) = zen;
}

void LidortRtDriver::setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1>& atmosphere_bb)
{
  Lidort_Fixed_Optical& foptical_inputs = lidort_interface_->lidort_fixin().optical();

  foptical_inputs.ts_surface_bb_input(surface_bb);

  // Thermal black body atmosphere inputs will be on levels instead of layers
  Range rlev(0, atmosphere_bb.extent(firstDim) - 1);
  Array<double, 1> thermal_bb_input( foptical_inputs.ts_thermal_bb_input() );
  thermal_bb_input(rlev) = atmosphere_bb;
}

void LidortRtDriver::setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                          const blitz::Array<double, 1>& ssa,
                                          const blitz::Array<double, 2>& pf)
{

  // Ranges for copying inputs to method
  Range rlay(0, od.extent(firstDim) - 1);
  Range rmom(0, pf.extent(firstDim) - 1);

  // Convienence references
  Lidort_Fixed_Optical& foptical_inputs = lidort_interface_->lidort_fixin().optical();
  Lidort_Modified_Optical& moptical_inputs = lidort_interface_->lidort_modin().moptical();

  // Vertical optical depth thicness values for all layers and threads
  Array<double, 1> deltau( foptical_inputs.ts_deltau_vert_input() );
  deltau(rlay) = od;

  // Single scattering albedos for all layers and threads
  Array<double, 1> omega( moptical_inputs.ts_omega_total_input() );
  omega(rlay) = where(ssa > 0.999, 0.999999, ssa);

  // For all layers n and threads t, Legrenre moments of
  // the phase function expansion multiplied by (2L+1);
  // initial value (L=0) should always be 1
  // phasmoms_total_input(n, L, t)
  // n = moments, L = layers, t = threads
  Array<double, 2> phasmoms( foptical_inputs.ts_phasmoms_total_input() );
  phasmoms(rmom, rlay) = where(abs(pf) > 1e-11, pf, 1e-11);
}

void LidortRtDriver::clear_linear_inputs()
{
  Lidort_Modified_Lincontrol& lincontrol = lidort_interface_->lidort_linmodin().mcont();

  // Flag for output of profile Jacobians
  lincontrol.ts_do_profile_linearization(false);

  // Flag for output of surface Jacobians
  lincontrol.ts_do_surface_linearization(false);
}

void LidortRtDriver::setup_linear_inputs(const ArrayAd<double, 1>& od,
                                         const ArrayAd<double, 1>& ssa,
                                         const ArrayAd<double, 2>& pf,
                                         bool do_surface_linearization)
{
  
  if(od.number_variable() > rt_pars_->max_atmoswfs()) {
    Exception err;
    err << "LIDORT has been compiled to allow a maximum of " << rt_pars_->max_atmoswfs()
        << " atmosphere derivatives to be calculated. We are trying to calculate "
        << od.number_variable() << " atmosphere derivatives";
    throw err;
  }

  Lidort_Fixed_Lincontrol& lincontrol = lidort_interface_->lidort_linfixin().cont();
  Lidort_Modified_Lincontrol& mlincontrol = lidort_interface_->lidort_linmodin().mcont();
  Lidort_Fixed_Linoptical& linoptical = lidort_interface_->lidort_linfixin().optical();

  // Number of profile weighting functions in layer n
  int natm_jac = od.number_variable();
  Array<int, 1> layer_jac_number( lincontrol.ts_layer_vary_number() );
  layer_jac_number = natm_jac;

  // Ranges for copying inputs to method
  Range ra(Range::all());
  Range rlay(0, od.rows() - 1);         // number of layers
  Range rjac(0, natm_jac - 1);

  Range rmom(0, pf.rows() - 1); // number phase function moments
  Range all(Range::all());

  // Flag for output of profile Jacobians
  // Unlikely we won't need these
  mlincontrol.ts_do_profile_linearization(true);

  // Flag for output of surface Jacobians
  // Certainly wouldn't need this if not retrieving ground
  mlincontrol.ts_do_surface_linearization(do_surface_linearization);

  // Enable surface black body jacobian when thermal emission is enabled
  mlincontrol.ts_do_surface_lbbf(do_thermal_emission);

  // Enable atmosphere black body jacobian when thermal emission is enabled
  mlincontrol.ts_do_atmos_lbbf(do_thermal_emission);

  // Check that we fit within the LIDORT configuration
  if(linoptical.ts_l_deltau_vert_input().cols() < od.rows()) {
    Exception e;
    e << "The number of layers you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << od.rows() << "\nand the maximum allowed is "
      << linoptical.ts_l_deltau_vert_input().cols() << "\n"
      << "\n"
      << "You might try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxlayer=value set to a larger value.\n";
    throw e;
  }
  if(linoptical.ts_l_deltau_vert_input().rows() < natm_jac) {
    Exception e;
    e << "The number of jacobians you are using exceeds the maximum allowed by\n"
      << "the current build of Lidort. The number requested is "
      << natm_jac << "\nand the maximum allowed is "
      << linoptical.ts_l_deltau_vert_input().rows() << "\n"
      << "\n"
      << "This number of jacobians is a function of the number of aerosols\n"
      << "in your state vector, so you can reduce the number of aerosols\n"
      << "\n"
      << "You might also try rebuilding with a larger value given to the configure\n"
      << "option --with-lidort-maxatmoswfs=value set to a larger value.\n";
    throw e;
  }

  // Setup optical linear inputs
  // LIDORT expects the following:
  // l_deltau = xi/tau * dtau/dxi
  // l_omega = xi/omega * domega/dxi
  // ... etc
  // The driver is handed dtau/dxi, domega/dxi ... etc
  // Where it is normally set up for xi being taug, taur, taua[0], taua[1]...
  //
  // LIDORT returns to us:
  // xi * dI/dxi
  // Therefore you will notice that by not multiplying dtau/dxi by xi and only
  // dividing by xi, we are cancelling out the xi in the result and hence
  // the driver really return dI/dxi
  Array<double, 2> l_deltau( linoptical.ts_l_deltau_vert_input()(rjac,rlay) );
  Array<double, 2> l_omega( linoptical.ts_l_omega_total_input()(rjac,rlay) );
  Array<double, 3> l_phasmoms( linoptical.ts_l_phasmoms_total_input()(rjac,rmom,rlay) );

  // Transpose these to match dimensions used internally
  l_deltau.transposeSelf(secondDim, firstDim);
  l_omega.transposeSelf(secondDim, firstDim);

  // Note that LIDORT does *not* take l_deltau etc. Rather it takes
  // what it calls "normalized derivative inputs". So for an optical
  // quantity like taug wrt to tau, this would be taug / tau *
  // dtau/dtaug.
  //
  // It then returns a normalized jacobian, which for taug would be
  // taug * d I /dtaug.
  //
  // Note that we actually leave out one of these factors, because it
  // effectively cancels out. We pass in 1 / tau * dtau / dtaug and
  // get back d I / dtaug.

  firstIndex i1; secondIndex i2;
  l_deltau = where(od.value()(i1) != 0, od.jacobian() / od.value()(i1), 0.0);

  Array<double, 1> ssa_limit(ssa.rows());
  ssa_limit = where(ssa.value() > 0.999, 0.999999, ssa.value());
  l_omega  = where(ssa_limit(i1) != 0, ssa.jacobian() / ssa_limit(i1), 0.0);

  if(pf.is_constant())
    l_phasmoms(rjac, rmom, rlay) = 0.0;
  else {
    blitz::Array<double, 2> pf_in( where(abs(pf.value()) > 1e-11, pf.value(), 1e-11) );

    // We need this loop since l_phasmoms and pf variables have jacobian data in different dimensions
    for (int jidx = 0; jidx < pf.number_variable(); jidx++)
      l_phasmoms(jidx, rmom, rlay) = pf.jacobian()(rmom, rlay, jidx) / pf_in(rmom, rlay)(i1,i2);
  }
}

/// Copy outputs from BRDF supplement into LIDORT Sup inputs types
void LidortRtDriver::copy_brdf_sup_outputs() const {

  // Copy BRDF outputs to LIDORT's BRDF inputs
  // Use LIDORT specific interfaces to ensure copying is done in memory correctly
  // Cannot mix VLIDORT and LIDORT BRDF objects here
  Brdf_Sup_Outputs& brdf_outputs = lidort_brdf_interface()->brdf_sup_out();
  Brdf_Linsup_Outputs& brdf_lin_outputs = lidort_brdf_interface()->brdf_linsup_out();

  Lidort_Sup_Brdf& lid_brdf = lidort_interface_->lidort_sup().brdf();
  Lidort_Linsup_Brdf& lid_lin_brdf = lidort_interface_->lidort_linsup().brdf();

  lid_brdf.copy_from_sup( brdf_outputs );
  lid_lin_brdf.copy_from_sup( brdf_lin_outputs );

}

void LidortRtDriver::calculate_rt() const
{
  // Must copy current BRDF supplement outputs into datastructures used by LIDORT
  copy_brdf_sup_outputs();

  // Call LIDORT for calculations
  lidort_interface_->run(false);
}

double LidortRtDriver::get_intensity() const
{
  // Total Intensity I(t,v,d,T) at output level t, output geometry v,
  // direction d
  return lidort_interface_->lidort_out().main().ts_intensity()(0,0, rt_pars_->upidx()-1);
}

void LidortRtDriver::copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_param, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const
{
  Lidort_Linatmos& lpoutputs = lidort_interface_->lidort_linout().atmos();
  Lidort_Linsurf& lsoutputs = lidort_interface_->lidort_linout().surf();

  Range ra(Range::all());

  // Surface Jacobians KR(r,t,v,d) with respect to surface variable r
  // at output level t, geometry v, direction d
  jac_surf_param.reference( lsoutputs.ts_surfacewf()(ra, 0, 0, rt_pars_->upidx()-1).copy() );

  // Get profile jacobians
  // Jacobians K(q,n,t,v,d) with respect to profile atmospheric variable
  // q in layer n, at output level t, geometry v, direction d
  jac_atm.reference( lpoutputs.ts_profilewf()(ra, ra, 0, 0, rt_pars_->upidx()-1).copy() );

  // Get surface temp jacobian if thermal emission is enabled
  if(do_thermal_emission) {
      jac_surf_temp = lsoutputs.ts_sbbwfs_jacobians()(0, 0, rt_pars_->upidx()-1);

      jac_atm_temp.reference( lpoutputs.ts_abbwfs_jacobians()(0, 0, ra, rt_pars_->upidx()-1).copy() );
  }
}
