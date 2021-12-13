#include "lidort_fixture.h"
#include "stokes_coefficient_constant.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

LidortDriverCommonFixture::LidortDriverCommonFixture() : sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 0.001;
  zen = 0.0;
  azm = 0.0;
  pure_nadir = false;
  do_multiple_scattering_only = false;

  nstreams = 4;
  nmoms = 2*nstreams;

  nlayer = 1;
  heights.resize(nlayer+1);
  od.resize(nlayer, 1);
  ssa.resize(nlayer, 1);
  pf.resize(nmoms, nlayer, 1);

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  Range all = Range::all();
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5; 

  surface_params = 0.0;
  surface_type = -1;
  od = 0.0;
  ssa = 0.0;

  taug = 0.0;
  taur = 0.0;
}

LidortDriverLambertianFixture::LidortDriverLambertianFixture() : LidortDriverCommonFixture()
{
  surface_params.resize(1); 
  surface_type = LAMBERTIAN;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir));  

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();
}

LidortDriverLambertianThermalFixture::LidortDriverLambertianThermalFixture() : LidortDriverCommonFixture()
{
  surface_params.resize(1); 
  surface_type = LAMBERTIAN;

  bool do_solar = false;
  bool do_thermal = true;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir, do_solar, do_thermal));

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();
}

LidortDriverCoxmunkFixture::LidortDriverCoxmunkFixture() : LidortDriverCommonFixture()
{
  surface_params.resize(4); 
  surface_type = COXMUNK;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir));

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();
}

LidortDriverBrdfVegFixture::LidortDriverBrdfVegFixture() : LidortDriverCommonFixture()
{
  surface_params.resize(5); 
  surface_type = BREONVEG;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir));

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();
}



/********************************************************************/

LidortRtCommonFixture::LidortRtCommonFixture() : sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 74.128288268999995;
  zen = 1.0e-6;
  azm = 12.504928589999992;
  pure_nadir = false;

  wn_arr.resize(2);
  wn_arr = 12929.94, 12930.30; 

  blitz::Array<double, 2> stokes_coef_v(3, 3);
  stokes_coef_v = 
    1,0,0,
    1,0,0,
    1,0,0;
  stokes_coefs.reset(new StokesCoefficientConstant(stokes_coef_v));
}

LidortLambertianFixture::LidortLambertianFixture(const std::string& Config_file) : 
  LidortRtCommonFixture(), LuaConfigurationFixture(Config_file)
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;
  
  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));
  
  pp_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_pp_plus_sscorrection";
}

LidortCoxmunkFixture::LidortCoxmunkFixture(const std::string& Config_file) : 
  LidortRtCommonFixture(), ConfigurationCoxmunkFixture(Config_file)
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));

  pp_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_pp_plus_sscorrection";

}

LidortCoxmunkPlusLambertianFixture::LidortCoxmunkPlusLambertianFixture() : LidortRtCommonFixture()
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));

  pp_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_pp_plus_sscorrection";

}

/*******************************************************************/

LidortLowHighCommon::LidortLowHighCommon() 
: sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 74.128288268999995;
  zen = 0;
  azm = 12.504928589999992;
  pure_nadir = true;

  blitz::Array<double, 2> stokes_coef_v(3, 3);
  stokes_coef_v = 
    1,0,0,
    1,0,0,
    1,0,0;
  stokes_coefs.reset(new StokesCoefficientConstant(stokes_coef_v));

  nstreams_low = 2;
  nmoms_low = 2*nstreams_low;

  nstreams_high = 8;
  nmoms_high = Lidort_Pars::instance().maxmoments_input;

}

LidortLowHighLambertianFixture::LidortLowHighLambertianFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));

}

LidortLowHighFixtureNoPolarization::LidortLowHighFixtureNoPolarization() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = false;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}

LidortLowHighCoxmunkFixture::LidortLowHighCoxmunkFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}

LidortLowHighCoxmunkPlusLambertianFixture::LidortLowHighCoxmunkPlusLambertianFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}
