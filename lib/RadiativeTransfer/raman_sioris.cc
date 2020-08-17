#include "raman_sioris.h"
#include "fp_serialize_support.h"

#include "unit.h"

#include "angle_util.h"

// for to_fortran
#include "linear_algebra.h"

using namespace blitz;
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RamanSiorisEffect::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectrumEffectImpBase)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure)
    & FP_NVP_(channel_index) & FP_NVP_(albedo) & FP_NVP_(padding_fraction)
    & FP_NVP_(do_upwelling) & FP_NVP_(jac_perturbation) & FP_NVP_(solar_zenith)
    & FP_NVP_(obs_zenith) & FP_NVP_(relative_azimuth) & FP_NVP_(scattering_angle)
    & FP_NVP_(temperature_layers) & FP_NVP_(atmosphere)
    & FP_NVP_(solar_model) & FP_NVP_(absorber);
}

FP_IMPLEMENT(RamanSiorisEffect);
#endif


extern "C" {
    void get_raman(int *nz, int *nw, int* maxnu, double *sza, double *vza, double *sca, double *albedo, bool *do_upwelling, const double *ts, const double *rhos, const double *wave, const double *sol, const double *taus, double *rspec, bool *problems);
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

    // Compute the maximum extent of internal get_raman arrays so we never run into
    // hitting a maximum value. The Fortran code computes an internal wavelenght grid
    // spaced every 1 wavenumber, but still in nm. This value represents the extent
    // it would expect.
    int maxnu = int(1e7/nm_grid_f(0)) - int(1e7/nm_grid_f(num_points-1)) + 2;

    Array<double, 1> raman_spec(num_points, ColumnMajorArray<1>());

    bool problems = false;

    get_raman(&num_layers, &num_points, &maxnu, &solar_zenith, &viewing_zenith, &scattering_angle, &albedo, &do_upwelling, temperature_layers.dataFirst(), air_number_density.dataFirst(), nm_grid_f.dataFirst(), solar_irradiance_f.dataFirst(), total_optical_depth_f.dataFirst(), raman_spec.dataFirst(), &problems);

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

RamanSiorisEffect::RamanSiorisEffect(double scale_factor,
                                     int channel_index,
                                     const DoubleWithUnit& solar_zenith, 
                                     const DoubleWithUnit& observation_zenith, 
                                     const DoubleWithUnit& relative_azimuth,
                                     const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
                                     const boost::shared_ptr<SolarModel>& solar_model,
                                     double albedo,
                                     double padding_fraction,
                                     bool do_upwelling,
                                     double jac_perturbation)
: SpectrumEffectImpBase(scale_factor),
  channel_index_(channel_index),
  albedo_(albedo),
  padding_fraction_(padding_fraction),
  do_upwelling_(do_upwelling),
  jac_perturbation_(jac_perturbation),
  atmosphere_(atmosphere),
  solar_model_(solar_model)
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
/// Extends the computed range out additional points to account for
/// the cut-off region in the Fortran model.
///-----------------------------------------------------------------------

void RamanSiorisEffect::apply_effect
(Spectrum& Spec,
 const ForwardModelSpectralGrid& UNUSED(Forward_model_grid)
) const
{
    Range ra = Range::all();

    // Convert dry air number density to the necessary units of molecules/cm^2
    Array<double, 1> dry_air_density = absorber_->dry_air_molecular_density_layer().value.value();
    dry_air_density *= 1e-4;

    // Compute a padded grid due to requirements of the fortran code
    // Pad 10% of the size of the input grid
    double pad_amount = 
        padding_fraction_ * (Spec.spectral_domain().data()(Spec.spectral_domain().rows()-1) - Spec.spectral_domain().data()(0));
    SpectralDomain padded_grid = Spec.spectral_domain().add_padding(DoubleWithUnit(pad_amount, Spec.spectral_domain().units()));

    // Compute total optical depth, hopefully use any caching that atmosphere class provides
    Array<double, 1> wn_grid = padded_grid.wavenumber();
    Array<double, 2> total_optical_depth(wn_grid.rows(), temperature_layers_.rows());

    for(int wn_idx = 0; wn_idx < wn_grid.rows(); wn_idx++) {
        total_optical_depth(wn_idx, ra) = atmosphere_->optical_depth_wrt_rt(wn_grid(wn_idx), channel_index_).value();
    }
    
    // Absolute magnitude of the solar spectrum does not seem to matter, no need to convert to a certain unit
    Array<double, 1> solar_spectrum = solar_model_->solar_spectrum(padded_grid).spectral_range().data();

    // Compute raman spectrum
    Array<double, 1> raman_spec(padded_grid.data().rows());
    raman_spec = compute_raman_sioris(solar_zenith_, obs_zenith_, scattering_angle_, albedo_, do_upwelling_, temperature_layers_, dry_air_density, padded_grid, solar_spectrum, total_optical_depth);

    // Load scale factor
    AutoDerivative<double> scale_factor = coefficient()(0);

    // Apply raman scattering scaling and compute jacobian through a perturbation
    ArrayAd<double, 1> spec_rad = Spec.spectral_range().data_ad();

    // Create spectrum with raman effect applied, index into padded raman scattering only for the points to be returned
    Array<double, 1> raman_applied_rad_val(spec_rad.rows());
    Array<double, 1> raman_applied_rad_pert(spec_rad.rows());

    for(int pad_idx = 0; pad_idx < padded_grid.data().rows(); pad_idx++) {
        int out_idx = padded_grid.sample_index()(pad_idx);
        if (out_idx >= 0) {
            raman_applied_rad_val(out_idx) = spec_rad.value()(out_idx) * (raman_spec(pad_idx) * scale_factor.value() + 1.0);
            raman_applied_rad_pert(out_idx) = spec_rad.value()(out_idx) * (raman_spec(pad_idx) * (scale_factor.value() + jac_perturbation_) + 1.0);
        }
    }

    // Save applied effect
    spec_rad.value() = raman_applied_rad_val;

    // Compute FD jacobian values
    if(!spec_rad.is_constant() && !coefficient().is_constant()) {
        for(int sv_idx = 0; sv_idx < spec_rad.jacobian().cols(); sv_idx++) {
            if(coefficient().jacobian()(0, sv_idx) > 0) {
                spec_rad.jacobian()(ra, sv_idx) = coefficient().jacobian()(0, sv_idx) * (raman_applied_rad_pert - raman_applied_rad_val) / jac_perturbation_;
            }
        }
    }

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
            new RamanSiorisEffect(coefficient().value()(0),
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
     << "  Solar zenith angle:     " << solar_zenith_ << "\n"
     << "  Observation angle:      " << obs_zenith_ << "\n"
     << "  Relative azimuth angle: " << relative_azimuth_ << "\n"
     << "  Scattering angle:       " << scattering_angle_ << "\n"
     << "  Albedo:                 " << albedo_ << "\n"
     << "  Do upwelling:           " << do_upwelling_ << "\n";
}
