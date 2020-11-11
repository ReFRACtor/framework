#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "generic_object_map.h"

#include "atmosphere_standard.h"
#include "atmosphere_legacy.h"
#include "altitude.h"
#include "fp_serialize_support.h"

#include <blitz/array.h>

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(atmosphere_standard, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IfstreamCs expected_data(test_data_dir() + "expected/atmosphere_standard/rt_parameters_each_layer");
  // Expected values were gotten by running the old Fortran code and
  // extracting out the answer from that.
  Array<double, 1> od_expect, ssa_expect;
  // Pick one level to check, because the matrix is big. No
  // significance in the level we selected here, it was just one in
  // the middle.
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_rt(12929.94, 0).value(), od_expect, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (atm->single_scattering_albedo_wrt_rt(12929.94, 0).value(), ssa_expect, 1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "O2").value(), 0.00018824861591009646, 1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "H2O").value(), 3.459804510441169e-06, 1e-4);
  Array<double, 2> scat_momsub
    (atm->phase_function_moments_wrt_rt(12929.94, 0, 4, 1).value()
     (Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub, 
                               scat_momsub_expect, 1e-4);
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_rt(12930.30, 0).value(), od_expect, 2e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_rt(12930.30, 0).value(), ssa_expect, 1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "O2").value(), 0.00019080438913330467, 1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "H2O").value(), 5.6828474797889968e-06, 1e-4);
  Array<double, 2> scat_momsub2(atm->phase_function_moments_wrt_rt(12930.30, 0, 4, 1).value()
                                (Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub2, 
                               scat_momsub_expect, 1e-4);
}

BOOST_AUTO_TEST_CASE(rayleigh_atmosphere)
{
  // There are three ways to get a Rayleigh atmosphere.
  // We check that leaving out the aerosol (null pointer),
  // defining a aerosol object with no aerosols defined, and
  // simply having the aerosol extinction set to 0 all give
  // the same results.

  boost::shared_ptr<AtmosphereStandard> atm_zeroext = atm->clone();
  // Set state vector so that the extinction coefficient of
  // atm_zeroext is 0. 
  StateVector sv;
  sv.add_observer(*atm_zeroext);
  atm_zeroext->attach_children_to_sv(sv);
  sv.update_state(config_initial_guess->initial_guess());
  blitz::Array<double, 1> sv_value = sv.state();
  for(int i = 0; i < sv.state_vector_name().rows(); ++i) {
    std::cerr << sv.state_vector_name()(i) << std::endl;
    if(sv.state_vector_name()(i).find("Aerosol Ext") != std::string::npos)
      sv_value(i) = 1e-20;
  }
  sv.update_state(sv_value);

  // Leave out the aerosol (null pointer)
  boost::shared_ptr<Pressure> pressure_clone = atm->pressure_ptr()->clone();
  boost::shared_ptr<Temperature> temperature_clone =
    atm->temperature_ptr()->clone();
  boost::shared_ptr<Rayleigh> rayleigh_clone = atm->rayleigh_ptr()->clone();
  boost::shared_ptr<Ground> ground_clone = atm->ground()->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr())
    alt_clone.push_back(a->clone());
  boost::shared_ptr<Absorber> absorber_clone =
    atm->absorber_ptr()->clone();
  boost::shared_ptr<RelativeHumidity> rh_clone =
    atm->relative_humidity_ptr()->clone();
  boost::shared_ptr<AerosolOptical> aerosol_null;
  boost::shared_ptr<AtmosphereStandard> 
    atm_rayleigh(new AtmosphereStandard(absorber_clone,
                                   pressure_clone,
                                   temperature_clone,
                                   rayleigh_clone,
                                   aerosol_null,
                                   rh_clone,
                                   ground_clone,
                                   alt_clone,
                                   atm->constant_ptr()));

  // Create an aerosol object with no aerosols defined
  std::vector<boost::shared_ptr<AerosolExtinction> > empty_aext;
  std::vector<boost::shared_ptr<AerosolProperty> >   empty_aprop;
  pressure_clone = atm->pressure_ptr()->clone();
  temperature_clone = atm->temperature_ptr()->clone();
  rayleigh_clone = atm->rayleigh_ptr()->clone();
  ground_clone = atm->ground()->clone();
  alt_clone.clear();
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr()) {
    alt_clone.push_back(a->clone());
  }
  absorber_clone =  atm->absorber_ptr()->clone();
  rh_clone = atm->relative_humidity_ptr()->clone();

  boost::shared_ptr<AerosolOptical> no_aerosol = boost::make_shared<AerosolOptical>(empty_aext,
                                                  empty_aprop,
                                                  pressure_clone,
                                                  rh_clone);
  boost::shared_ptr<AtmosphereStandard>
    atm_no_aerosol(new AtmosphereStandard(absorber_clone,
                                   pressure_clone,
                                   temperature_clone,
                                   rayleigh_clone,
                                   no_aerosol,
                                   rh_clone,
                                   ground_clone,
                                   alt_clone,
                                   atm->constant_ptr()));

  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->ground()->surface_parameter(12929.94, 0).value(),
                           atm_zeroext->ground()->surface_parameter(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->ground()->surface_parameter(12929.94, 0).value(),
                           atm_no_aerosol->ground()->surface_parameter(12929.94, 0).value());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_spectrometer(),
                    atm_zeroext->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_spectrometer(),
                    atm_no_aerosol->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_layer(),
                    atm_zeroext->number_layer());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_layer(),
                    atm_no_aerosol->number_layer());

  std::ofstream of_rayleigh("atm_rayleigh.new");
  of_rayleigh << *atm_rayleigh << std::endl;

  std::ofstream of_zeroext("atm_zeroext.new");
  of_zeroext << *atm_zeroext << std::endl;

  std::cerr << "atm_rayleigh = " << atm_rayleigh->optical_depth_wrt_rt(12929.94, 0).value() << std::endl
            << "atm_zeroext = " << atm_zeroext->optical_depth_wrt_rt(12929.94, 0).value() << std::endl;
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->optical_depth_wrt_rt(12929.94, 0).value(),
                           atm_zeroext->optical_depth_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->optical_depth_wrt_rt(12929.94, 0).value(),
                           atm_no_aerosol->optical_depth_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->single_scattering_albedo_wrt_rt(12929.94, 0).value(), 
                           atm_zeroext->single_scattering_albedo_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->single_scattering_albedo_wrt_rt(12929.94, 0).value(),
                           atm_no_aerosol->single_scattering_albedo_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "O2").value(), 
     atm_zeroext->column_optical_depth(12929.94, 0, "O2").value(), 
     1e-6);
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "O2").value(),
     atm_no_aerosol->column_optical_depth(12929.94, 0, "O2").value(),
     1e-6);
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "H2O").value(),
     atm_zeroext->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "H2O").value(),
     atm_no_aerosol->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);
  Range r1(0,atm_rayleigh->phase_function_moments_wrt_rt(12738.853381475927, 0).rows() - 1);
  Range r2(0,
           atm_rayleigh->phase_function_moments_wrt_rt(12738.853381475927, 0).depth() - 1);
  Range ra = Range::all();
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->phase_function_moments_wrt_rt(12929.94, 0).value(),
                           atm_zeroext->phase_function_moments_wrt_rt(12929.94, 0).value()(r1, ra, r2));
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->phase_function_moments_wrt_rt(12929.94, 0).value(),
                           atm_no_aerosol->phase_function_moments_wrt_rt(12929.94, 0).value()(r1, ra, r2));
}

BOOST_AUTO_TEST_CASE(uplooking_atmosphere)
{
  // We check that leaving out the ground gives the same results 
  boost::shared_ptr<Pressure> pressure_clone = atm->pressure_ptr()->clone();
  boost::shared_ptr<Temperature> temperature_clone =
    atm->temperature_ptr()->clone();
  boost::shared_ptr<Rayleigh> rayleigh_clone = atm->rayleigh_ptr()->clone();
  boost::shared_ptr<Aerosol> aerosol_clone =
    atm->aerosol_ptr()->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr())
    alt_clone.push_back(a->clone());
  boost::shared_ptr<Absorber> absorber_clone =
    atm->absorber_ptr()->clone();
  boost::shared_ptr<RelativeHumidity> rh_clone =
    atm->relative_humidity_ptr()->clone();
  boost::shared_ptr<Ground> ground_null;
  boost::shared_ptr<AtmosphereStandard> 
    atm_uplooking(new AtmosphereStandard(absorber_clone,
                                    pressure_clone,
                                    temperature_clone,
                                    rayleigh_clone,
                                    aerosol_clone,
                                    rh_clone,
                                    ground_null,
                                    alt_clone,
                                    atm->constant_ptr()));
  boost::shared_ptr<StateVector> sv_uplooking(new StateVector());
  attach_atmosphere_to_sv(atm_uplooking, sv_uplooking);

  BOOST_CHECK_EQUAL(atm_uplooking->number_spectrometer(),
                    atm->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_uplooking->number_layer(),
                    atm->number_layer());

  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->optical_depth_wrt_rt(12929.94, 0).value(),
                           atm->optical_depth_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->single_scattering_albedo_wrt_rt(12929.94, 0).value(),
                           atm->single_scattering_albedo_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_CLOSE
    (atm_uplooking->column_optical_depth(12929.94, 0, "O2").value(), 
     atm->column_optical_depth(12929.94, 0, "O2").value(), 1e-6);
  BOOST_CHECK_CLOSE
    (atm_uplooking->column_optical_depth(12929.94, 0, "H2O").value(),
     atm->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);

  ArrayAd<double, 3> pf_up(atm_uplooking->phase_function_moments_wrt_rt(12929.94, 0));
  ArrayAd<double, 3> pf_gnd(atm->phase_function_moments_wrt_rt(12929.94, 0));
  BOOST_CHECK_MATRIX_CLOSE(pf_up.value(), pf_gnd.value());
  BOOST_CHECK_MATRIX_CLOSE(pf_up.jacobian(), pf_gnd.jacobian());

  Array<double, 2> l_opdel_uplooking = atm_uplooking->optical_depth_wrt_rt(12929.94, 0).jacobian();
  Array<double, 2> l_opdel = atm->optical_depth_wrt_rt(12929.94, 0).jacobian();

  Array<double, 2> l_ssa_uplooking = atm_uplooking->single_scattering_albedo_wrt_rt(12929.94, 0).jacobian();
  Array<double, 2> l_ssa = atm->single_scattering_albedo_wrt_rt(12929.94, 0).jacobian();

  BOOST_CHECK_MATRIX_CLOSE(l_opdel_uplooking, l_opdel);
  BOOST_CHECK_MATRIX_CLOSE(l_ssa_uplooking, l_ssa);
}

BOOST_AUTO_TEST_CASE(no_aerosols)
{
  // We check that leaving out the aerosol gives the same results as
  // having no aerosols particles defined.


  // Create an AerosolOptical object with no aerosols in it
  boost::shared_ptr<Pressure> pressure_clone = atm->pressure_ptr()->clone();
  boost::shared_ptr<RelativeHumidity> rh_clone = atm->relative_humidity_ptr()->clone();
  std::vector<boost::shared_ptr<AerosolExtinction> > empty_aext;
  std::vector<boost::shared_ptr<AerosolProperty> >   empty_aprop;
  boost::shared_ptr<AerosolOptical> no_aerosol = boost::make_shared<AerosolOptical>(empty_aext,
                                                  empty_aprop,
                                                  pressure_clone,
                                                  rh_clone);

  boost::shared_ptr<Temperature> temperature_clone =
    atm->temperature_ptr()->clone();
  boost::shared_ptr<Rayleigh> rayleigh_clone = atm->rayleigh_ptr()->clone();
  boost::shared_ptr<Ground> ground_clone = atm->ground()->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr())
    alt_clone.push_back(a->clone());
  boost::shared_ptr<Absorber> absorber_clone =
    atm->absorber_ptr()->clone();

  boost::shared_ptr<AtmosphereStandard>
    atm_no_aerosol(new AtmosphereStandard(absorber_clone,
                                   pressure_clone,
                                   temperature_clone,
                                   rayleigh_clone,
                                   no_aerosol,
                                   rh_clone,
                                   ground_clone,
                                   alt_clone,
                                   atm->constant_ptr()));

  boost::shared_ptr<AerosolOptical> aerosol_null;
  boost::shared_ptr<AtmosphereStandard>
    atm_rayleigh(new AtmosphereStandard(absorber_clone,
                                   pressure_clone,
                                   temperature_clone,
                                   rayleigh_clone,
                                   aerosol_null,
                                   rh_clone,
                                   ground_clone,
                                   alt_clone,
                                   atm->constant_ptr()));
  /*

  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->ground()->surface_parameter(12929.94, 0).value(),
                           no_aerosol->ground()->surface_parameter(12929.94, 0).value());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_spectrometer(),
                    no_aerosol->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_layer(),
                    no_aerosol->number_layer());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->optical_depth_wrt_rt(12929.94, 0).value(),
                           no_aerosol->optical_depth_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->single_scattering_albedo_wrt_rt(12929.94, 0).value(),
                           no_aerosol->single_scattering_albedo_wrt_rt(12929.94, 0).value());
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "O2").value(),
     no_aerosol->column_optical_depth(12929.94, 0, "O2").value(),
     1e-6);
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "H2O").value(),
     no_aerosol->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);
  Range r1(0,atm_rayleigh->phase_function_moments_wrt_rt(12738.853381475927, 0).rows() - 1);
  Range r2(0,
           atm_rayleigh->phase_function_moments_wrt_rt(12738.853381475927, 0).depth() - 1);
  Range ra = Range::all();
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->phase_function_moments_wrt_rt(12929.94, 0).value(),
                           no_aerosol->phase_function_moments_wrt_rt(12929.94, 0).value()(r1, ra, r2));
                           */
}

BOOST_AUTO_TEST_CASE(legacy)
{
    double test_wn = 13179.0;
    int test_chan = 0;

    boost::shared_ptr<AtmosphereLegacy> atm_legacy
        (new AtmosphereLegacy(atm->absorber_ptr(), atm->pressure_ptr(), atm->temperature_ptr(), atm->rayleigh_ptr(),
                              atm->aerosol_ptr(), atm->relative_humidity_ptr(), atm->ground(), 
                              atm->altitude_ptr(), atm->constant_ptr()));

    // Total optical depth
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_state_vector(test_wn, test_chan).value(), 
                                 atm->optical_depth_wrt_state_vector(test_wn, test_chan).value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_state_vector(test_wn, test_chan).jacobian(), 
                                 atm->optical_depth_wrt_state_vector(test_wn, test_chan).jacobian(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_rt(test_wn, test_chan).value(), 
                                 atm->optical_depth_wrt_rt(test_wn, test_chan).value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->optical_depth_wrt_rt(test_wn, test_chan).jacobian(), 
                                 atm->optical_depth_wrt_rt(test_wn, test_chan).jacobian(), 1e-10);

    // Column optical depth
    BOOST_CHECK_CLOSE(atm_legacy->column_optical_depth(test_wn, test_chan, "O2").value(), 
                      atm->column_optical_depth(test_wn, test_chan, "O2").value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->column_optical_depth(test_wn, test_chan, "O2").gradient(), 
                                 atm->column_optical_depth(test_wn, test_chan, "O2").gradient(), 1e-10);

    ArrayAd<double, 2> gas_od_per_part(atm->optical_properties(test_wn, test_chan)->gas_optical_depth_per_particle());

    // Single scattering albedo
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).value(), 
                                 atm->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).jacobian(),
                                 atm->single_scattering_albedo_wrt_state_vector(test_wn, test_chan).jacobian(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_rt(test_wn, test_chan).value(), 
                                 atm->single_scattering_albedo_wrt_rt(test_wn, test_chan).value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(atm_legacy->single_scattering_albedo_wrt_rt(test_wn, test_chan).jacobian(),
                                 atm->single_scattering_albedo_wrt_rt(test_wn, test_chan).jacobian(), 1e-10);

    // Total pf which include rayleigh portion which is not exposed by Atmosphere interface
    ArrayAd<double, 3> tot_pf_wrt_sv_expt_1 = atm_legacy->phase_function_moments_wrt_state_vector(test_wn, test_chan);
    ArrayAd<double, 3> tot_pf_wrt_sv_calc_1 = atm->phase_function_moments_wrt_state_vector(test_wn, test_chan);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_sv_expt_1.value(), tot_pf_wrt_sv_calc_1.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_sv_expt_1.jacobian(), tot_pf_wrt_sv_calc_1.jacobian(), 1e-10);

    ArrayAd<double, 3> tot_pf_wrt_rt_expt_1 = atm_legacy->phase_function_moments_wrt_rt(test_wn, test_chan);
    ArrayAd<double, 3> tot_pf_wrt_rt_calc_1 = atm->phase_function_moments_wrt_rt(test_wn, test_chan);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_rt_expt_1.value(), tot_pf_wrt_rt_calc_1.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_rt_expt_1.jacobian(), tot_pf_wrt_rt_calc_1.jacobian(), 1e-10);

    // Check that phase function can be recomputed with a different number of moments and scattering
    ArrayAd<double, 3> tot_pf_wrt_sv_expt_2 = atm_legacy->phase_function_moments_wrt_state_vector(test_wn, test_chan, 200, 1);
    ArrayAd<double, 3> tot_pf_wrt_sv_calc_2 = atm->phase_function_moments_wrt_state_vector(test_wn, test_chan, 200, 1);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_sv_expt_2.value(), tot_pf_wrt_sv_calc_2.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_sv_expt_2.jacobian(), tot_pf_wrt_sv_calc_2.jacobian(), 1e-10);

    ArrayAd<double, 3> tot_pf_wrt_rt_expt_2 = atm_legacy->phase_function_moments_wrt_rt(test_wn, test_chan, 200, 1);
    ArrayAd<double, 3> tot_pf_wrt_rt_calc_2 = atm->phase_function_moments_wrt_rt(test_wn, test_chan, 200, 1);

    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_rt_expt_2.value(), tot_pf_wrt_rt_calc_2.value(), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(tot_pf_wrt_rt_expt_2.jacobian(), tot_pf_wrt_rt_calc_2.jacobian(), 1e-10);

}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
  boost::shared_ptr<GenericObjectMap> m =
    boost::make_shared<GenericObjectMap>();
  (*m)["atm"] = atm;
  (*m)["statev" ] = statev;
  std::string d = serialize_write_string(m);
  if(false)
    std::cerr << d;
  boost::shared_ptr<GenericObjectMap> mr =
    serialize_read_string<GenericObjectMap>(d);
  boost::shared_ptr<AtmosphereStandard> atmr =
    mr->get<AtmosphereStandard>("atm");
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_rt(12929.94, 0).value(),
			       atmr->optical_depth_wrt_rt(12929.94, 0).value(),
			       1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (atm->single_scattering_albedo_wrt_rt(12929.94, 0).value(),
     atmr->single_scattering_albedo_wrt_rt(12929.94, 0).value(),
     1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "O2").value(),
		    atmr->column_optical_depth(12929.94, 0, "O2").value(),
		    1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "H2O").value(),
		    atmr->column_optical_depth(12929.94, 0, "H2O").value(),
		    1e-4);
  Array<double, 2> scat_momsub
    (atm->phase_function_moments_wrt_rt(12929.94, 0, 4, 1).value()
     (Range::all(), 9, Range::all()));
  Array<double, 2> scat_momsub2
    (atmr->phase_function_moments_wrt_rt(12929.94, 0, 4, 1).value()
     (Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub, scat_momsub2, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_rt(12930.30, 0).value(),
			       atmr->optical_depth_wrt_rt(12930.30, 0).value(),
			       2e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (atm->single_scattering_albedo_wrt_rt(12930.30, 0).value(),
     atmr->single_scattering_albedo_wrt_rt(12930.30, 0).value(), 1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "O2").value(),
		    atmr->column_optical_depth(12930.30, 0, "O2").value(),
		    1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "H2O").value(),
		    atmr->column_optical_depth(12930.30, 0, "H2O").value(),
		    1e-4);
  boost::shared_ptr<StateVector> statevr =
    mr->get<StateVector>("statev");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(atmosphere_standard_jac, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  IfstreamCs expected_data(test_data_dir() + 
                           "expected/atmosphere_standard/rt_parameters_each_layer");
  Array<double, 1> od_expect, ssa_expect;
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(
           atm.optical_depth_wrt_state_vector(12929.94, 0).value(), 
           od_expect, 1e-6);
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 1> od = 
    atm.optical_depth_wrt_state_vector(wn, spec_index);
  Array<double, 1> od0(od.shape());
  od0 = od.value();
  Array<double, 2> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(od0.shape());
    jacfd = (atm.optical_depth_wrt_rt(wn,spec_index).value() - od0) 
      / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(ssa_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  IfstreamCs expected_data(test_data_dir() + 
                           "expected/atmosphere_standard/rt_parameters_each_layer");
  Array<double, 1> od_expect, ssa_expect;
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(
    atm.single_scattering_albedo_wrt_state_vector(12929.94, 0).value(), 
    ssa_expect, 1e-6);
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 1> ssa = 
    atm.single_scattering_albedo_wrt_state_vector(wn, spec_index);
  Array<double, 1> ssa0(ssa.shape());
  ssa0 = ssa.value();
  Array<double, 2> jac = ssa.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(ssa0.shape());
    jacfd = (atm.single_scattering_albedo_wrt_rt(wn,spec_index).value() - ssa0) 
      / epsilon(i);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
        std::cerr << i << ": " << jac(Range::all(), i) << "\n"
                  << jacfd << "\n"
                  << diff << "\n";
    }
    // There are a wide range in the size of the difference, because
    // the Jacobian values vary wildly in size from one row to the
    // next. We looked at all the values, and things look
    // reasonable. To have an automated test, we look at having either
    // the difference being very small, or the relative difference
    // small.
    Array<double, 1> diff(jac(Range::all(), i) - jacfd);
    BOOST_CHECK(max(abs(diff)) < 4e-6 || max(abs(where(jacfd == 0, 0, diff / jacfd))) < 5e-4);
  }
}

BOOST_AUTO_TEST_CASE(phase_function_moments_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 3> sm = 
    atm.phase_function_moments_wrt_state_vector(wn, spec_index);
  Array<double, 3> sm0(sm.shape());
  sm0 = sm.value();
  Array<double, 4> jac = sm.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 3> jacfd(sm.shape());
    jacfd = (atm.phase_function_moments_wrt_rt(wn,spec_index).value() - sm0) 
      / epsilon(i);
    Array<double, 3> diff(jac(Range::all(), Range::all(), Range::all(), i) - 
                          jacfd);
    if(false) {                        // Can turn this off to dump values,
                                // if needed for debugging
      if(max(abs(diff)) > 0) {
        std::cerr << i << ": " << max(abs(diff)) << " " 
                  << max(abs(where(abs(jacfd) < 1e-15, 0, diff / jacfd))) 
                  << "\n";
        std::cerr << jac(Range(0,9), Range::all(), 0, i) << "\n"
                  << jacfd(Range(0,9), Range::all(), 0) << "\n";
      }
    }
  BOOST_CHECK(all(abs(diff) < 4e-6 || abs(where(jacfd == 0, 0, diff / jacfd)) < 1e-4));
  }
}

BOOST_AUTO_TEST_CASE(optical_depth_timing)
{
  is_timing_test();
  RtAtmosphere& atm = *config_atmosphere;
  int i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01) {
    ArrayAd<double, 1> od = atm.optical_depth_wrt_rt(wn, 0);
    if(++i % 1000 == 0)
      std::cerr << "Done with " << i << "\n"
                << atm.timer_info() << "\n";
  }
}

BOOST_AUTO_TEST_SUITE_END()
