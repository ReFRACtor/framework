#include "unit.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"

#include "hdf_file.h"
#include "raman_sioris.h"

#include "standard_forward_model.h"
#include "simple_fixed_spectrum_sampling.h"
#include "ground_lambertian.h"
#include "default_constant.h"
#include "absorber_absco.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(raman_sioris, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(offline_data)
{
    // Unit test data was converted from test_raman_input.dat file from MUSES OMI code in idl-retrieve
    // that was originally provided by Xiong Liu
    HdfFile offline_data(test_data_dir() + "/in/raman_sioris/offline_test_data.h5");

    int UNUSED(nw) = offline_data.read_field<int>("/nw");
    int UNUSED(nz) = offline_data.read_field<int>("/nz");
    double sza = offline_data.read_field<double>("/sza");
    double vza = offline_data.read_field<double>("/vza");
    double sca = offline_data.read_field<double>("/sca");

    Array<double, 1> ts = offline_data.read_field<double, 1>("/ts");
    Array<double, 1> rhos = offline_data.read_field<double, 1>("/rhos");

    Array<double, 1> wave = offline_data.read_field<double, 1>("/wave");
    Array<double, 1> sol = offline_data.read_field<double, 1>("/sol");

    Array<double, 2> taus = offline_data.read_field<double, 2>("/taus");

    bool do_upwelling = true;
    double albedo = 0.5;

    SpectralDomain wave_sd(wave, units::nm);

    Array<double, 1> rspec_calc = FullPhysics::compute_raman_sioris(sza, vza, sca, albedo, do_upwelling, ts, rhos, wave_sd, wave_sd, sol, taus);

    Array<double, 1> rspec_expt = offline_data.read_field<double, 1>("rspec");

    BOOST_CHECK_MATRIX_CLOSE_TOL(rspec_expt, rspec_calc, 1e-6);

}

BOOST_AUTO_TEST_CASE(effect)
{
    int channel_idx = 0;
    double scale_factor = 1.9;

    // Convert this so we can access ->spectrum_effect
    boost::shared_ptr<StandardForwardModel> fm = boost::dynamic_pointer_cast<StandardForwardModel>(config_forward_model);

    // This is highly dependent on the solar model being the SpectrumEffect at this index
    boost::shared_ptr<SolarModel> solar_model = boost::dynamic_pointer_cast<SolarModel>(fm->spectrum_effect()[0][0]);

    // Cast to the object we need
    boost::shared_ptr<AtmosphereStandard> atm = boost::dynamic_pointer_cast<AtmosphereStandard>(config_atmosphere);

    // Construct a rayleigh only atmosphere since we only need rayleigh scattering for Raman scattering
    boost::shared_ptr<Constant> constant(new DefaultConstant());
    boost::shared_ptr<AbsorberAbsco> absorber(
            new AbsorberAbsco(std::vector<boost::shared_ptr<AbsorberVmr> >(),
                              atm->pressure_ptr(), atm->temperature_ptr(), atm->altitude_ptr(), 
                              std::vector<boost::shared_ptr<GasAbsorption> >(), constant));

    boost::shared_ptr<AtmosphereStandard> atm_rayleigh(
            new AtmosphereStandard(absorber, atm->pressure_ptr(), atm->temperature_ptr(), atm->rayleigh_ptr(),
                                   atm->relative_humidity_ptr(), atm->ground(), atm->altitude_ptr(), constant));
    auto g = boost::dynamic_pointer_cast<GroundLambertian>(atm->ground());
    // Set albedo to 1.0, to match expected test results.
    g->sub_state_vector_values()(0) = 1.0;
    g->sub_state_vector_values()(1) = 0;
    g->sub_state_vector_values()(2) = 1.0;
    g->sub_state_vector_values()(3) = 0;
    g->sub_state_vector_values()(4) = 1.0;
    g->sub_state_vector_values()(5) = 0;

    // Pull out angles we need
    DoubleWithUnit solar_zenith = config_level_1b->solar_zenith(channel_idx);
    DoubleWithUnit observation_zenith = config_level_1b->sounding_zenith(channel_idx);
    DoubleWithUnit relative_azimuth = config_level_1b->relative_azimuth(channel_idx);

    // Make a grid from 740 to 770 nm
    // The grid must be wide enough for the RamanSioris code to calculate any values
    Array<double, 1> grid_vals(30);
    double wave_nm = 740;
    for(int grid_idx = 0; grid_idx < grid_vals.rows(); grid_idx++) {
        grid_vals(grid_idx) = wave_nm++;
    }
    SpectralDomain sd = SpectralDomain(grid_vals, units::nm);
    double pad_amount =
      0.1 * (sd.data()(sd.rows()-1) - sd.data()(0));
    SpectralDomain padded_grid = sd.add_padding(DoubleWithUnit(pad_amount, units::nm));
    RamanSiorisEffect raman =
      RamanSiorisEffect(padded_grid, scale_factor, channel_idx, solar_zenith,
			observation_zenith, relative_azimuth,
			atm_rayleigh, solar_model);

    ForwardModelSpectralGrid
      fm_spec_grid(config_instrument, config_spectral_window,
		   config_spectrum_sampling);

    int num_jac = 10;
    ArrayAd<double, 1> spec_range(sd.data().rows(), num_jac);
    spec_range.value() = 1.0;
    spec_range.jacobian() = 0.0;

    // Units don't matter here, but lets just assign something reasonable
    Unit rad_units("ph / s / m^2 / micron W / (cm^-1) / (ph / (s) / (micron)) sr^-1");
    Spectrum spec(sd, SpectralRange(spec_range, rad_units));

    // Set up statevector stuff so that we can properly
    // test the jacobians going through our created fluorescence
    // class
    StateVector sv;
    sv.add_observer(raman);
    Array<double,1> x(num_jac);
    x(Range(0,0)) = scale_factor;
    sv.update_state(x);

    // Apply effect
    raman.apply_effect(spec, fm_spec_grid);

    Array<double, 1> applied_effect = spec.spectral_range().data();
    Array<double, 2> calc_jac( spec.spectral_range().data_ad().jacobian() );

    // These are not special vals, just values computed after other testing looked okay and captured
    // for automatic testing. In fact this set up is probably not optimal as this feature should
    // be wider and not have so much null values on either end
    Array<double, 1> expt_vals(30); 
    expt_vals =
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1.0409, 1.07519, 1.15625, 1.1122, 1.07699, 1.06549, 1.02736, 1.01291, 1.00533, 1.00195, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

    BOOST_CHECK_MATRIX_CLOSE_TOL(expt_vals, applied_effect, 1e-5);

    // Check jacobian by doing and outside perturbation and comparing
    double perturbation = 0.001;

    ArrayAd<double, 1> spec_range_pert(sd.data().rows(), 0);
    spec_range_pert.value() = 1.0;
    Spectrum spec_pert(sd, SpectralRange(spec_range_pert, spec.spectral_range().units()));

    x(Range(0,0)) = scale_factor + perturbation;
    sv.update_state(x);

    raman.apply_effect(spec_pert, fm_spec_grid);

    Array<double, 1> pert_effect = spec_pert.spectral_range().data();

    Array<double, 1> expt_jac( (pert_effect - applied_effect) / perturbation );
    BOOST_CHECK_MATRIX_CLOSE_TOL(expt_jac, calc_jac(Range::all(), 0), 1e-10);

}

BOOST_AUTO_TEST_SUITE_END()
