#include "vlidort_fixture.h"
#include "stokes_coefficient_constant.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

VLidortDriverCommonFixture::VLidortDriverCommonFixture(int Nstokes)
  : nstokes(Nstokes), sza(3), zen(3), azm(3)
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
  // Third dimension carries the 6 de-Rooij scattering matrix elements
  // [a1, a2, a3, a4, b1, b2]
  pf.resize(nmoms, nlayer, 6, 1);

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // Rayleigh scattering phase matrix moments with depolarization factor = 0.
  // These are the de-Rooij expansion coefficients used by both
  // RayleighGreekMoment and the VLIDORT tester (PROBLEM_RAY).
  Range all = Range::all();
  pf = 0.0;
  double depol = 0.0;
  double beta2 = (1.0 - depol) / (2.0 + depol);
  pf(0, all, 0) = 1.0;                                  // a1, L=0
  pf(2, all, 0) = beta2;                                // a1, L=2
  pf(2, all, 1) = 6.0 * beta2;                          // a2, L=2
  pf(1, all, 3) = 3.0 * (1.0 - 2.0*depol) / (2.0 + depol); // a4, L=1
  pf(2, all, 4) = std::sqrt(6.0) * beta2;               // b1, L=2

  surface_params = 0.0;
  surface_type = -1;
  od = 0.0;
  ssa = 0.0;

  taug = 0.0;
  taur = 0.0;
}

VLidortDriverLambertianFixture::VLidortDriverLambertianFixture(int Nstokes)
  : VLidortDriverCommonFixture(Nstokes)
{
  surface_params.resize(1);
  surface_type = LAMBERTIAN;

  vlidort_driver.reset(new VLidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir, nstokes));

  // Turn off delta-m scaling
  vlidort_driver->vlidort_interface()->vlidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  vlidort_driver->set_plane_parallel();
}

VLidortDriverLambertianThermalFixture::VLidortDriverLambertianThermalFixture(int Nstokes)
  : VLidortDriverCommonFixture(Nstokes)
{
  surface_params.resize(1);
  surface_type = LAMBERTIAN;

  bool do_solar = false;
  bool do_thermal = true;

  vlidort_driver.reset(new VLidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir, nstokes, do_solar, do_thermal));

  // Turn off delta-m scaling
  vlidort_driver->vlidort_interface()->vlidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  vlidort_driver->set_plane_parallel();
}

VLidortDriverCoxmunkFixture::VLidortDriverCoxmunkFixture(int Nstokes)
  : VLidortDriverCommonFixture(Nstokes)
{
  surface_params.resize(4);
  surface_type = COXMUNK;

  vlidort_driver.reset(new VLidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir, nstokes));

  // Turn off delta-m scaling
  vlidort_driver->vlidort_interface()->vlidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  vlidort_driver->set_plane_parallel();
}

VLidortDriverVectorFixture::VLidortDriverVectorFixture(int Nstokes)
  : VLidortDriverCommonFixture(Nstokes)
{
  // Off-nadir, non-zero azimuth geometry so Q and U are non-trivial.
  // Must match test_vlidort_refractor.f90 in the offline tester.
  sza = 45.0;
  zen = 35.0;
  azm = 60.0;

  surface_params.resize(1);
  surface_type = LAMBERTIAN;

  vlidort_driver.reset(new VLidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir, nstokes));

  // Turn off delta-m scaling
  vlidort_driver->vlidort_interface()->vlidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  vlidort_driver->set_plane_parallel();
}
