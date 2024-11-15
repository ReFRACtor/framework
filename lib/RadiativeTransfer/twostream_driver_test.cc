#include "twostream_driver.h"
#include "unit_test_support.h"
#include "ground.h"
#include "old_constant.h"
#include "lidort_driver.h"
#include "spurr_brdf_types.h"
#include "planck.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(twostream_driver, GlobalFixture)

void test_twostream(boost::shared_ptr<TwostreamRtDriver>& twostream_driver, ArrayAd<double, 1>& surface_params, ArrayAd<double, 1>& taug, ArrayAd<double, 1>& taur, ArrayAd<double, 1>& taua, Array<double, 1>& pert_atm, Array<double, 1>& pert_surf, bool do_solar, bool do_thermal, bool debug_output)
{
  blitz::Array<double, 1> sza, zen, azm;
  sza.resize(1); zen.resize(1); azm.resize(1);
  sza = 30.0;
  zen = 20.0;
  azm = 0.0;
  bool pure_nadir = false;
  int nmoms = 2;

  // Aerosol properties
  double aer_prop_ssa = 0.9;
  double aer_prop_asym = 0.7;
  double depol = 0.0;

  int nlayer = taug.rows();
  int nparam = taug.number_variable();

  // Set up convenience ranges
  Range all = Range::all();
  Range rjac(0,nparam-1);
  Range rlay(0,nlayer-1);
  Range rsurf(0,surface_params.number_variable()-1);
  double refl_ts;
  double refl_lid;

  blitz::Array<double, 2> jac_atm_ts;
  blitz::Array<double, 1> jac_surf_param_ts;
  double jac_surf_temp_ts;
  blitz::Array<double, 1> jac_atm_temp_ts;

  blitz::Array<double, 2> jac_atm_lid;
  blitz::Array<double, 1> jac_surf_param_lid;
  double jac_surf_temp_lid;
  blitz::Array<double, 1> jac_atm_temp_lid;

  // Use LIDORT for comparison
  int lid_nstreams = 1;
  int lid_nmoms = 2*lid_nstreams;
  bool do_multiple_scattering_only = true;
  LidortRtDriver lidort_driver =
    LidortRtDriver(lid_nstreams, lid_nmoms, 
                   do_multiple_scattering_only,
                   twostream_driver->surface_type(), zen, pure_nadir,
                   do_solar, do_thermal);  

  // Plane-parallel
  lidort_driver.set_plane_parallel();

  // Simple height grid evenly spaced
  Array<double, 1> heights(nlayer+1);
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // Set up atmosphere jacobians to be one for taug and the other for taur
  taug.jacobian()(all, 0) = 1;
  taur.jacobian()(all, 1) = 1;

  if (pert_atm.rows() > 2 && pert_atm(2) != 0.0) {
      taua.jacobian()(all, 2) = 1;
  }

  ArrayAd<double, 1> od(nlayer, nparam);
  ArrayAd<double, 1> ssa(nlayer, nparam);
  
  ArrayAd<double,1> ray_wt(nlayer, nparam);
  ArrayAd<double,1> aer_wt(nlayer, nparam);

  for(int lay_idx = 0; lay_idx < nlayer; lay_idx++) {
    od(lay_idx) = taur(lay_idx) + taug(lay_idx) + taua(lay_idx);
    ssa(lay_idx) = (taur(lay_idx) + aer_prop_ssa*taua(lay_idx)) / od(lay_idx);

    ray_wt(lay_idx) = taur(lay_idx) / (taur(lay_idx) + aer_prop_ssa * taua(lay_idx));
    aer_wt(lay_idx) = 1.0 - ray_wt(lay_idx);
  } 

  // Compute phase function
  ArrayAd<double, 2> pf(nmoms+1, nlayer, nparam);
  pf = 0.0;

  for(int lay_idx = 0; lay_idx < nlayer; lay_idx++) {
      pf(0, lay_idx) = 1.0;
      pf(2, lay_idx) = ray_wt(lay_idx) * ( (1.0 - depol) / (2.0 - depol) );

      for(int mom_idx = 1; mom_idx <= nmoms; mom_idx++) {
          pf(mom_idx, lay_idx) = pf(mom_idx, lay_idx) + aer_wt(lay_idx) * (2*mom_idx+1) * pow(aer_prop_asym, mom_idx);
      }
  }

  // Set up thermal inputs if enabled
  double bb_surface = 0.0;
  Array<double, 1> bb_atm(nlayer + 1);
  if (do_thermal) {
      // Just some nominal values to calculate the planck function
      double wn = 568.69;
      double temperature = 290.0;
      bb_surface = planck(wn, temperature);
      bb_atm = planck(wn, temperature);
  }

  // Set up jacobians
  surface_params.jacobian() = 1.0;

  // Make copies of surface params as it will be modified
  // internally by the routines
  ArrayAd<double, 1> lid_surface_params;
  lid_surface_params = surface_params;
  ArrayAd<double, 1> ts_surface_params;
  ts_surface_params = surface_params;

  // Run lidort and 2stream to generate values for comparison
  lidort_driver.reflectance_and_jacobian_calculate
    (heights, sza(0), zen(0), azm(0),
     twostream_driver->surface_type(), lid_surface_params,
     od, ssa, pf, refl_lid, 
     jac_atm_lid, jac_surf_param_lid, jac_surf_temp_lid, jac_atm_temp_lid,
     bb_surface, bb_atm);
  twostream_driver->reflectance_and_jacobian_calculate
    (heights, sza(0), zen(0), azm(0),
     twostream_driver->surface_type(), ts_surface_params,
     od, ssa, pf, refl_ts, 
     jac_atm_ts, jac_surf_param_ts, jac_surf_temp_ts, jac_atm_temp_ts,
     bb_surface, bb_atm);

  if(debug_output) {
    std::cerr << "refl_lid = " << refl_lid << std::endl
              << "refl_ts  = " << refl_ts << std::endl;
  }

  BOOST_CHECK_CLOSE(refl_lid, refl_ts, 7e-5);

  blitz::Array<double, 2> jac_atm_fd(nparam, nlayer);
  jac_atm_fd = 0.0;

  double refl_fd;

  for(int l_idx = 0; l_idx < nlayer; l_idx++) {
    for(int p_idx = 0; p_idx < pert_atm.rows(); p_idx++) {
      blitz::Array<double,1> taug_pert( taug.value().copy() );
      blitz::Array<double,1> taur_pert( taur.value().copy() );
      blitz::Array<double,1> taua_pert( taua.value().copy() );

      blitz::Array<double,1> od_pert( od.value().shape() );
      blitz::Array<double,1> ssa_pert( ssa.value().shape() );
      blitz::Array<double,2> pf_pert( pf.value().copy() );
   
     switch (p_idx) {
      case 0:
        taug_pert(l_idx) += pert_atm(p_idx);
        break;
      case 1:
        taur_pert(l_idx) += pert_atm(p_idx);
        break;
      case 2:
        taua_pert(l_idx) += pert_atm(p_idx);
        break;
      }

      od_pert = taur_pert + taug_pert + taua_pert;
      ssa_pert = (taur_pert + aer_prop_ssa*taua_pert) / od_pert;

      double ray_wt_pert = taur_pert(l_idx) / (taur_pert(l_idx) + aer_prop_ssa * taua_pert(l_idx));
      double aer_wt_pert = 1.0 - ray_wt_pert;

      pf_pert(all, l_idx) = 0;
      pf_pert(0, l_idx) = 1.0;
      pf_pert(2, l_idx) = ray_wt_pert * ( (1.0 - depol) / (2.0 - depol) );

      for(int mom_idx = 1; mom_idx <= nmoms; mom_idx++) {
          pf_pert(mom_idx, l_idx) = pf_pert(mom_idx, l_idx) + aer_wt_pert * (2*mom_idx+1) * pow(aer_prop_asym, mom_idx);
      }

      refl_fd = twostream_driver->reflectance_calculate
        (heights, sza(0), zen(0), azm(0),
         twostream_driver->surface_type(), surface_params.value().copy(),
         od_pert, ssa_pert, pf_pert,
         bb_surface, bb_atm);

      if (pert_atm(p_idx) != 0.0) {
        jac_atm_fd(p_idx, l_idx) = (refl_fd - refl_ts) / pert_atm(p_idx);
      }
    }
  }

  if(debug_output) {
    std::cerr << setprecision(8)
              << "jac_atm_lid = " << jac_atm_lid(rjac,rlay).transpose(1,0)
              << "jac_atm_fd = " << jac_atm_fd(rjac, rlay).transpose(1,0)
              << "jac_atm_ts = " << jac_atm_ts(rjac,rlay).transpose(1,0) << std::endl;
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_lid(rjac,rlay), jac_atm_ts(rjac,rlay), 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_ts(rjac,rlay), jac_atm_fd(rjac,rlay), 5e-4);

  if(blitz::any(surface_params.value() > 0.0)) {
    // Check surface jacobians against finite difference
    blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param_ts.rows() );
    jac_surf_param_fd = 0.0;

    for(int p_idx = 0; p_idx < pert_surf.rows(); p_idx++) {
      blitz::Array<double,1> surface_params_pert( surface_params.rows() );
      surface_params_pert = surface_params.value();
      surface_params_pert(p_idx) += pert_surf(p_idx);

      refl_fd = twostream_driver->reflectance_calculate
        (heights, sza(0), zen(0), azm(0),
         twostream_driver->surface_type(), surface_params_pert,
         od.value(), ssa.value(), pf.value(),
         bb_surface, bb_atm);

      jac_surf_param_fd(p_idx) = (refl_fd - refl_ts) / pert_surf(p_idx);

      // Adjust analytic jacobians to have same meaning as finite difference one
      jac_surf_param_lid(p_idx) *= lid_surface_params.jacobian()(p_idx, 0);
      jac_surf_param_ts(p_idx) *= ts_surface_params.jacobian()(p_idx, 0);
    }

    if(debug_output) {
      std::cerr << "jac_surf_param_lid = " << jac_surf_param_lid(rsurf) << std::endl
                << "jac_surf_param_fd = " << jac_surf_param_fd << std::endl
                << "jac_surf_param_ts = " << jac_surf_param_ts << std::endl;
    }

    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param_lid(rsurf), jac_surf_param_ts(rsurf), 1e-7);
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param_ts, jac_surf_param_fd, 1e-7);
  }
}

void test_twostream_lambertian(boost::shared_ptr<TwostreamRtDriver>& twostream_driver, bool do_solar, bool do_thermal, bool debug_output)
{
  int nlayer = 2;
  int nparam = 3;
  ArrayAd<double, 1> surface_params(1, 1); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);
  ArrayAd<double, 1> taua(nlayer, nparam);

  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(nparam);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(1);

  ////////////////
  // Surface only
  if (do_thermal && !do_solar) {
      // Make emissivity a large number
      surface_params.value() = 1.0e-6;
  } else {
      surface_params.value() = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0;

  if (debug_output) std::cerr << "Surface only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4, 0;
  pert_surf = 1e-8;
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

  ////////////////
  // Rayleigh only

  if (do_thermal && !do_solar) {
      // For thermal mode make the surface makes emissivity close to 0. Finite difference goes screwy if this is set == 1.0
      surface_params.value() = 0.98;
  } else {
      // 2stream will divide by 0 in thermal mode if surface param is set to identically 0 for lambertian mode
      surface_params.value() = 1.0e-6;
  }

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0;

  if (debug_output) std::cerr << "Rayleigh only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-3, -1e-4, 0;
  pert_surf = 1e-8;
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

  ////////////////
  // Gas + Surface
  if (do_thermal && !do_solar) {
      surface_params.value() = 1.0e-6;
  } else {
      surface_params.value() = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;
  taua = 0;

  if (debug_output) std::cerr << "Gas + Surface" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4, 0;
  pert_surf = 1e-8;
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

  ////////////////////////////
  // Gas + Surface + Aerosol
  if (do_thermal && !do_solar) {
      surface_params.value() = 1.0e-6;
  } else {
      surface_params.value() = 1.0;
  }

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;
  taua = 1.0e-2/nlayer;

  if (debug_output) std::cerr << "Gas + Surface + Aerosol" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4, 1e-4;
  pert_surf = 1e-8;
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);

}

BOOST_AUTO_TEST_CASE(lambertian_solar)
{
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false;
  int surface_type = LAMBERTIAN;
  int nlayer = 2;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);


  test_twostream_lambertian(twostream_driver, do_solar, do_thermal,
                            debug_output);
}

BOOST_AUTO_TEST_CASE(lambertian_solar_serialization)
{
  if(!have_serialize_supported())
    return;
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false;
  int surface_type = LAMBERTIAN;
  int nlayer = 2;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);


  std::string d = serialize_write_string(twostream_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<TwostreamRtDriver> tdriver_r =
    serialize_read_string<TwostreamRtDriver>(d);

  test_twostream_lambertian(tdriver_r, do_solar, do_thermal,
                            debug_output);
}

BOOST_AUTO_TEST_CASE(lambertian_thermal)
{
  bool do_solar = false;
  bool do_thermal = true;
  bool debug_output = false;
  int surface_type = LAMBERTIAN;
  int nlayer = 2;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);

  test_twostream_lambertian(twostream_driver, do_solar, do_thermal,
                            debug_output);
}

BOOST_AUTO_TEST_CASE(lambertian_thermal_serialization)
{
  if(!have_serialize_supported())
    return;
  bool do_solar = false;
  bool do_thermal = true;
  bool debug_output = false;
  int surface_type = LAMBERTIAN;
  int nlayer = 2;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);

  std::string d = serialize_write_string(twostream_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<TwostreamRtDriver> tdriver_r =
    serialize_read_string<TwostreamRtDriver>(d);

  test_twostream_lambertian(tdriver_r, do_solar, do_thermal,
                            debug_output);
}

BOOST_AUTO_TEST_CASE(coxmunk)
{
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false; 

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(4, 3); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);
  ArrayAd<double, 1> taua(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = COXMUNK;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);
  pert_atm = 1e-4, 1e-4;
  
  // Surface only
  surface_params(0) = 7;
  surface_params(1) = 1.334;
  surface_params(2) = 0.8;
  surface_params(3) = 0.0;

  // Surface perturbations
  blitz::Array<double, 1> pert_surf(surface_params.number_variable());
  pert_surf = 1e-4, 1e-4, 1e-4;
  //pert_surf(1) = sqrt(surface_params(1).value()*surface_params(1).value() + pert_surf(1)) - surface_params(1).value();

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0.0;

  if (debug_output) std::cerr << "Coxmunk" << std::endl << "----------------------------" << std::endl;  
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_CASE(coxmunk_serialization)
{
  if(!have_serialize_supported())
    return;
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false; 

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(4, 3); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);
  ArrayAd<double, 1> taua(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = COXMUNK;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);
  pert_atm = 1e-4, 1e-4;
  
  // Surface only
  surface_params(0) = 7;
  surface_params(1) = 1.334;
  surface_params(2) = 0.8;
  surface_params(3) = 0.0;

  // Surface perturbations
  blitz::Array<double, 1> pert_surf(surface_params.number_variable());
  pert_surf = 1e-4, 1e-4, 1e-4;
  //pert_surf(1) = sqrt(surface_params(1).value()*surface_params(1).value() + pert_surf(1)) - surface_params(1).value();

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0.0;

  if (debug_output) std::cerr << "Coxmunk" << std::endl << "----------------------------" << std::endl;  
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);

  std::string d = serialize_write_string(twostream_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<TwostreamRtDriver> tdriver_r =
    serialize_read_string<TwostreamRtDriver>(d);
  
  test_twostream(tdriver_r, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_CASE(brdf)
{
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false;

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(5, 1); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);
  ArrayAd<double, 1> taua(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = BREONVEG;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(5);
  pert_surf = 1e-8;

  ////////////////
  // Surface only
  surface_params.value()(0) = 1.0; // rahman kernel factor
  surface_params.value()(1) = 0.1; // hotspot parameter
  surface_params.value()(2) = 0.3; // asymmetry
  surface_params.value()(3) = 1.5; // anisotropy_parameter
  surface_params.value()(4) = 1.0; // breon kernel factor

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0.0;

  pert_atm = 1e-4, 1e-4;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);
  test_twostream(twostream_driver, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_CASE(brdf_serialization)
{
  if(!have_serialize_supported())
    return;
  bool do_solar = true;
  bool do_thermal = false;
  bool debug_output = false;

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(5, 1); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);
  ArrayAd<double, 1> taua(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = BREONVEG;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(5);
  pert_surf = 1e-8;

  ////////////////
  // Surface only
  surface_params.value()(0) = 1.0; // rahman kernel factor
  surface_params.value()(1) = 0.1; // hotspot parameter
  surface_params.value()(2) = 0.3; // asymmetry
  surface_params.value()(3) = 1.5; // anisotropy_parameter
  surface_params.value()(4) = 1.0; // breon kernel factor

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;
  taua = 0.0;

  pert_atm = 1e-4, 1e-4;
  boost::shared_ptr<TwostreamRtDriver> twostream_driver =
    boost::make_shared<TwostreamRtDriver>(nlayer, surface_type, false,
                                          do_solar, do_thermal);

  // Plane-parallel
  twostream_driver->twostream_interface()->do_plane_parallel(true);

  std::string d = serialize_write_string(twostream_driver);
  if(false)
    std::cerr << d;
  boost::shared_ptr<TwostreamRtDriver> tdriver_r =
    serialize_read_string<TwostreamRtDriver>(d);
  
  test_twostream(tdriver_r, surface_params, taug, taur, taua, pert_atm, pert_surf, do_solar, do_thermal, debug_output);
}

BOOST_AUTO_TEST_CASE(valgrind_problem)
{
  // This catchs a particular valgrind problem we encountered in the
  // full physics code. This originally came from wavenumber
  // 13050.47 in the first iteration of tccon_sounding_1. There are
  // other valgrind errors there, but the hope is fixing this one will
  // fix all of them.

  // Note this has been fixed, but we'll leave this test in place
  // for now.
  
  IfstreamCs captured_input(test_data_dir() + "expected/twostream/valgrind_problem");
  int nlayer, nparm, surface_type, do_fullquadrature;
  captured_input >> nlayer >> nparm >> surface_type 
                 >> do_fullquadrature;
  TwostreamRtDriver d(nlayer, surface_type, (do_fullquadrature == 1));
                      
  double refl;
  Array<double, 2> jac_atm;
  Array<double, 1> jac_surf_param;
  double jac_surf_temp;
  blitz::Array<double, 1> jac_atm_temp;
  Array<double, 1> height;
  ArrayAd<double, 1> surface_param, od, ssa;
  ArrayAd<double, 2> pf;
  double sza, zen, azm;
  captured_input >> height >> sza >> azm >> zen >> surface_type
                 >> surface_param >> od >> ssa >> pf;
  d.reflectance_and_jacobian_calculate(height, sza, azm, zen, 
                                       surface_type, surface_param,
                                       od, ssa, pf, refl,
                                       jac_atm, jac_surf_param, jac_surf_temp, jac_atm_temp);
  // This will trigger an error when we run with valgrind. Note that
  // we *don't* actually see a NaN here, rather this conditional
  // triggers the unitialized value error
  for(int i = 0; i < od.jacobian().rows(); ++i)
    for(int j = 0; j < od.jacobian().cols(); ++j) {
      if(std::isnan(jac_atm(j,i)))
        std::cerr << "Nan at jac_atm(" << j << ", " << i << ")\n";
    }
  for(int i = 0; i < jac_surf_param.rows(); ++i)
    if(std::isnan(jac_surf_param(i)))
      std::cerr << "Nan at jac_surf_param(" << i << ")\n";
}
BOOST_AUTO_TEST_SUITE_END()
