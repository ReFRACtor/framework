#include "raman_sioris.h"

#include "unit.h"

#include "angle_util.h"

// for to_fortran
#include "linear_algebra.h"

using namespace blitz;
using namespace FullPhysics;

extern "C" {
    void get_raman(int *nz, int *nw, double *sza, double *vza, double *sca, double *albedo, bool *do_upwelling, const double *ts, const double *rhos, const double *wave, const double *sol, const double *taus, double *rspec, bool *problems);
}

//-----------------------------------------------------------------------
/// Interface wrapper around the Sioris Fortran Raman scattering code.
/// This wrapper ensures that dimensions between arrays match and converts
/// arrays to column first ordering.
//-----------------------------------------------------------------------

Array<double, 1> FullPhysics::compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const Array<double, 1> &temperature_layers, const Array<double, 1>& air_number_density, const SpectralDomain &grid, const Array<double, 1> &solar_irradiance, const Array<double, 2> &total_optical_depth)
{

    int num_layers = temperature_layers.rows(); // nz
    int num_points = grid.data().rows(); // nw

    if (air_number_density.rows() != num_layers) {
        Exception err;
        err << "Number of air_number_density layers: " << air_number_density.rows() 
            << " does not match the number of layers specified by the temperature layers: " << num_layers;
        throw err;
    }
    
    if (total_optical_depth.cols() != num_layers) {
        Exception err;
        err << "Number of total_optical_depth layers: " << total_optical_depth.rows() 
            << " does not match the number of layers specified by the temperature layers: " << num_layers;
        throw err;
    }

    if(solar_irradiance.rows() != num_points) {
        Exception err;
        err << "Number of solar_irradiance grid points: " << solar_irradiance.rows() 
            << " does not match the number of grid points: " << num_points;
        throw err;
    }

    if(total_optical_depth.rows() != num_points) {
        Exception err;
        err << "Number of total_optical_depth grid points: " << total_optical_depth.rows() 
            << " does not match the number of grid points: " << num_points;
        throw err;
    }

    Array<double, 1> nm_grid = grid.convert_wave(units::nm);

    // Check if we need to reverse the grid because we are coming from wavenumbers
    // Reverse all relevant arrays since the fortran code assumes increasing ordering
    // in the spectral domain
    bool reverse_grid = nm_grid(0) > nm_grid(nm_grid.rows()-1);

    Array<double, 1> nm_grid_f(num_points, ColumnMajorArray<1>());
    Array<double, 2> total_optical_depth_f(num_points, num_layers, ColumnMajorArray<2>());
    Array<double, 1> solar_irradiance_f(num_points, ColumnMajorArray<1>());

    if (reverse_grid) {
        // Make local arrays with a view of the same memory as passed in, then
        // reverse them which twiddles the indexing internally. We can not twiddle
        // the passed in arrays directly since they are const
        Array<double, 2> total_optical_depth_reversed(total_optical_depth);
        total_optical_depth_reversed.reverseSelf(firstDim);

        Array<double, 1> solar_irradiance_reversed(solar_irradiance);
        solar_irradiance_reversed.reverseSelf(firstDim);

        // Copy values into memory seen by fortan in their reverse order
        nm_grid_f = nm_grid.reverse(firstDim);
        total_optical_depth_f = total_optical_depth_reversed;
        solar_irradiance_f = solar_irradiance_reversed;
    } else {
        nm_grid_f = nm_grid;
        total_optical_depth_f = total_optical_depth;
        solar_irradiance_f = solar_irradiance;
    }

    Array<double, 1> raman_spec(num_points, ColumnMajorArray<1>());

    bool problems = false;

    get_raman(&num_layers, &num_points, &solar_zenith, &viewing_zenith, &scattering_angle, &albedo, &do_upwelling, temperature_layers.dataFirst(), air_number_density.dataFirst(), nm_grid_f.dataFirst(), solar_irradiance_f.dataFirst(), total_optical_depth_f.dataFirst(), raman_spec.dataFirst(), &problems);

    if (problems) {
        throw Exception("raman_sioris fortran code encountered an internal issue.");
    }

    // Reverse back to order of calling grid
    if (reverse_grid) {
        raman_spec.reverseSelf(firstDim);
    }

    return raman_spec;

}

//-----------------------------------------------------------------------
/// Create a new Raman application SpectrumEffect using the Sioris 
/// method:
/// Sioris, C. E., & Evans, W. F. (2000).
/// Impact of rotational Raman scattering in the O2A band. 
/// Geophysical research letters, 27(24), 4085-4088.
/// 
/// Requires angles of the observation. Uses atmosphere for computing
/// optical depths. Uses the solar model for solar irradiance values.
/// The solar irradiance values must be commensurate with the expected
/// units of Ph/s/cm^2/nn
/// 
/// The retrieved albedo is not used so as to not tie this model to a
/// certain surface model. The value of albedo does not seem to have
/// a significant enough effect that an approximate value can suffice.
//-----------------------------------------------------------------------

RamanSiorisEffect::RamanSiorisEffect(double scale_factor, bool used_flag,
                                     int channel_index,
                                     const DoubleWithUnit& solar_zenith, 
                                     const DoubleWithUnit& observation_zenith, 
                                     const DoubleWithUnit& relative_azimuth,
                                     const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
                                     const boost::shared_ptr<SolarModel>& solar_model,
                                     double albedo,
                                     bool do_upwelling)
: SpectrumEffectImpBase(scale_factor, used_flag),
  atmosphere_(atmosphere), solar_model_(solar_model), channel_index_(channel_index), albedo_(albedo), do_upwelling_(do_upwelling)
{
    // Convert angles to degrees since these should not change
    solar_zenith_ = solar_zenith.convert(units::deg).value;
    obs_zenith_ = observation_zenith.convert(units::deg).value;
    relative_azimuth_ = relative_azimuth.convert(units::deg).value; // stored for use in cloning
    scattering_angle_ = scattering_angle(observation_zenith, solar_zenith, relative_azimuth).convert(units::deg).value;

    absorber_ = boost::dynamic_pointer_cast<AbsorberAbsco>(atmosphere_->absorber_ptr());

    if(!absorber_) {
        throw Exception("Expected absorber type from AtmosphereStandard to be of type AbsorberAbsco");
    }

    // Initialize the temperature levels
    compute_temp_layers(*atmosphere_->pressure_ptr());

    // Attach ourselves to get updates from the Pressure class to recompute the temperature layers whenever
    // the pressure updates
    atmosphere_->pressure_ptr()->add_observer(*this);

}

//-----------------------------------------------------------------------
/// Applies the effect of Raman scattering to the passed spectrum.
/// The effect is multiplicative.
///
/// TODO: Extend the computed range out additional points to account for
/// the cut-off region in the fortran model. Until this is done,
/// a couple nm on either size will not have an effect applied.
//-----------------------------------------------------------------------

void RamanSiorisEffect::apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const
{
    Range ra = Range::all();

    // Convert dry air number density to the necessary units of molecules/cm^2
    Array<double, 1> dry_air_density = absorber_->dry_air_molecular_density_layer().value.value();
    dry_air_density *= 1e-4;

    // Compute total optical depth, hopefully use any caching that atmosphere class provides
    Array<double, 1> wn_grid = Spec.spectral_domain().wavenumber();
    Array<double, 2> total_optical_depth(wn_grid.rows(), temperature_layers_.rows());

    for(int wn_idx = 0; wn_idx < wn_grid.rows(); wn_idx++) {
        total_optical_depth(wn_idx, ra) = atmosphere_->optical_depth_wrt_rt(wn_grid(wn_idx), channel_index_).value();
    }
    
    // Convert to the expect solar spectrum units Ph / s / cm^2 / nm^1
    Array<double, 1> solar_spectrum = solar_model_->solar_spectrum(Spec.spectral_domain()).spectral_range().convert(Unit("Ph s^-1 cm^-2 nm^-1")).data();

    // Compute raman spectrum
    Array<double, 1> raman_spec(wn_grid.rows());
    raman_spec = compute_raman_sioris(solar_zenith_, obs_zenith_, scattering_angle_, albedo_, do_upwelling_, temperature_layers_, dry_air_density, Spec.spectral_domain(), solar_spectrum, total_optical_depth);

    // Load scale factor
    AutoDerivative<double> scale_factor = coefficient()(0);

    // Apply raman scattering scaling and compute jacobian through a perturbation
    ArrayAd<double, 1> spec_rad = Spec.spectral_range().data_ad();

    Array<double, 1> raman_applied_rad(spec_rad.value() * (raman_spec * scale_factor.value() + 1.0));
    //Array<double, 1> raman_applied_rad_pert = spec_rad.data() * (raman_spec * (scale_factor.value() + perturbation_) + 1.0);

    spec_rad.value() = raman_applied_rad;

}

//-----------------------------------------------------------------------
/// Computes the per layer temperature grid using a layer pressure grid
/// consisting of average pressure values per layer.
//-----------------------------------------------------------------------

void RamanSiorisEffect::compute_temp_layers(const Pressure& pressure)
{
    ArrayAd<double, 1> pgrid_lev = pressure.pressure_grid().value;

    temperature_layers_.resize(pgrid_lev.rows()-1);

    for(int lay_idx = 0; lay_idx < temperature_layers_.rows(); lay_idx++) {
        AutoDerivativeWithUnit<double> press_lay( (pgrid_lev(lay_idx) + pgrid_lev(lay_idx+1)) / 2, pressure.pressure_grid().units );

        temperature_layers_(lay_idx) = atmosphere_->temperature_ptr()->temperature(press_lay).value.value();
    }
}

//-----------------------------------------------------------------------
/// Create a new copy of this class using the same internal state 
/// as it currently exists
//-----------------------------------------------------------------------

boost::shared_ptr<SpectrumEffect> RamanSiorisEffect::clone() const
{

    return boost::shared_ptr<SpectrumEffect>( 
            new RamanSiorisEffect(coefficient().value()(0), used_flag_value()(0),
                                  channel_index_,
                                  DoubleWithUnit(solar_zenith_, units::deg),
                                  DoubleWithUnit(obs_zenith_, units::deg),
                                  DoubleWithUnit(relative_azimuth_, units::deg),
                                  atmosphere_,
                                  solar_model_,
                                  albedo_,
                                  do_upwelling_));
}

//-----------------------------------------------------------------------
/// Output a textual representation of the class
//-----------------------------------------------------------------------

void RamanSiorisEffect::print(std::ostream& Os) const
{
  Os << "RamanSiorisEffect" << std::endl
     << "  Scale Factor:           " << coefficient().value()(0) << "\n"
     << "  Retrieval flag:         " << used_flag_value()(0) << "\n"
     << "  Solar zenith angle:     " << solar_zenith_ << "\n"
     << "  Observation angle:      " << obs_zenith_ << "\n"
     << "  Relative azimuth angle: " << relative_azimuth_ << "\n"
     << "  Scattering angle:       " << scattering_angle_ << "\n"
     << "  Albedo:                 " << albedo_ << "\n"
     << "  Do upwelling:           " << do_upwelling_ << "\n";
}
