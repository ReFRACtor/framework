#ifndef VLIDORT_FIXTURE_H
#define VLIDORT_FIXTURE_H

#include "serialized_configuration_fixture.h"
#include "vlidort_rt.h"

namespace FullPhysics {

/****************************************************************//**
 Base fixture for testing VLidortRtDriver. Mirrors
 LidortDriverCommonFixture but carries the 6 de-Rooij scattering
 matrix elements in the third dimension of the phase function and a
 configurable number of Stokes parameters.
*******************************************************************/

class VLidortDriverCommonFixture {
public:
  VLidortDriverCommonFixture(int nstokes = 3);
  virtual ~VLidortDriverCommonFixture() {};

  /// Number of Stokes parameters
  int nstokes;

  /// Vlidort driver to use in testing
  boost::shared_ptr<VLidortRtDriver> vlidort_driver;

  /// solar zenith, view zenith, and azimuth angles
  blitz::Array<double, 1> sza, zen, azm;
  bool pure_nadir;
  bool do_multiple_scattering_only;

  int nstreams;
  int nmoms;

  int nlayer;
  blitz::Array<double, 1> heights;

  blitz::Array<double, 1> surface_params;
  int surface_type;

  ArrayAd<double, 1> od;
  ArrayAd<double, 1> ssa;
  ArrayAd<double, 3> pf;

  double taug;
  double taur;
};

/****************************************************************//**
 Fixture for VLidortRtDriver lambertian mode
*******************************************************************/

class VLidortDriverLambertianFixture : public VLidortDriverCommonFixture {
public:
  VLidortDriverLambertianFixture(int nstokes = 3);
  virtual ~VLidortDriverLambertianFixture() {};
};

/****************************************************************//**
 Fixture for VLidortRtDriver coxmunk mode
*******************************************************************/

class VLidortDriverCoxmunkFixture : public VLidortDriverCommonFixture {
public:
  VLidortDriverCoxmunkFixture(int nstokes = 3);
  virtual ~VLidortDriverCoxmunkFixture() {};
};

/****************************************************************//**
 Fixture for VLidortRtDriver in a polarizing off-nadir geometry with a
 Lambertian surface, used to validate the vector (nstokes = 3, 4) Stokes
 output against the offline VLIDORT reference driver
 (reference/rt/vlidort_2p8p3_lin_tester test_vlidort_refractor.f90).
*******************************************************************/

class VLidortDriverVectorFixture : public VLidortDriverCommonFixture {
public:
  VLidortDriverVectorFixture(int nstokes = 3);
  virtual ~VLidortDriverVectorFixture() {};
};

}
#endif
