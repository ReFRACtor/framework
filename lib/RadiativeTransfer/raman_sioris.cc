#include "raman_sioris.h"
#include "fp_serialize_support.h"

#include "unit.h"

#include "angle_util.h"

// for to_fortran
#include "linear_algebra.h"

using namespace blitz;
using namespace FullPhysics;

// The "edge" we need to the desired range of the Raman calculation
const double RamanSiorisEffect::raman_edge_wavenumber = 218;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void RamanSiorisEffect::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpectrumEffectImpBase)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObserverPressure)
    & FP_NVP_(solar_and_odepth_spec_domain) & FP_NVP_(channel_index)
    & FP_NVP_(do_upwelling) & FP_NVP_(solar_zenith)
    & FP_NVP_(obs_zenith) & FP_NVP_(relative_azimuth)
    & FP_NVP_(scattering_angle) & FP_NVP_(atmosphere)
    & FP_NVP_(solar_model);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void RamanSiorisEffect::save(Archive &UNUSED(ar),
			     const unsigned int UNUSED(version)) const
{
}

template<class Archive>
void RamanSiorisEffect::load(Archive &UNUSED(ar),
			     const unsigned int UNUSED(version))
{
  compute_temp_layers(*atmosphere_->pressure_ptr());
}

FP_IMPLEMENT(RamanSiorisEffect);
#endif


extern "C" {
  void get_raman(int *nz, int *nw, int *nw_out, int* maxnu, double *sza, double *vza, double *sca, double *albedo, bool *do_upwelling, const double *ts, const double *rhos, const double *wave, double *wave_out, const double *sol, const double *taus, double *rspec, bool *problems);
}

//-----------------------------------------------------------------------
/// Interface wrapper around the Sioris Fortran Raman scattering code.
/// This wrapper ensures that dimensions between arrays match and converts
/// arrays to column first ordering.
//-----------------------------------------------------------------------

Array<double, 1> FullPhysics::compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const Array<double, 1> &temperature_layers, const Array<double, 1>& air_number_density, const SpectralDomain &grid, const SpectralDomain &grid_out, const Array<double, 1> &solar_irradiance, const Array<double, 2> &total_optical_depth)
{

    int num_layers = temperature_layers.rows(); // nz
    int num_points = grid.data().rows(); // nw
    int num_points_out = grid_out.data().rows(); // nw_out

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
    Array<double, 1> nm_grid_out = grid_out.convert_wave(units::nm);

    // Check if we need to reverse the grid because we are coming from wavenumbers
    // Reverse all relevant arrays since the fortran code assumes increasing ordering
    // in the spectral domain
    bool reverse_grid = nm_grid(0) > nm_grid(nm_grid.rows()-1);
    bool reverse_grid_out = nm_grid_out(0) > nm_grid_out(nm_grid_out.rows()-1);

    Array<double, 1> nm_grid_f(num_points, ColumnMajorArray<1>());
    Array<double, 1> nm_grid_out_f(num_points_out, ColumnMajorArray<1>());
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

    if (reverse_grid_out) {
        // Copy values into memory seen by fortan in their reverse order
        nm_grid_out_f = nm_grid_out.reverse(firstDim);
    } else {
        nm_grid_out_f = nm_grid_out;
    }
    
    // Compute the maximum extent of internal get_raman arrays so we never run into
    // hitting a maximum value. The Fortran code computes an internal wavelenght grid
    // spaced every 1 wavenumber, but still in nm. This value represents the extent
    // it would expect.
    int maxnu = int(1e7/nm_grid_f(0)) - int(1e7/nm_grid_f(num_points-1)) + 2;

    // Also check that we cover at least the minimum grid size. If we
    // don't the Fortran code will seg fault, and we don't have enough
    // data to calculate anything anyways.
    if((maxnu - 2) < RamanSiorisEffect::raman_edge_wavenumber) {
      Exception e;
      e << "The Solar_and_odepth_spec_domain is too small for the raman code\n"
	<< "  Width in wavenumber: " << maxnu - 2 << "\n"
	<< "  Minimum wavenumber required: "
	<< RamanSiorisEffect::raman_edge_wavenumber << "\n";
      throw e;
    }
    
    Array<double, 1> raman_spec(num_points_out, ColumnMajorArray<1>());

    bool problems = false;

    get_raman(&num_layers, &num_points, &num_points_out, &maxnu, &solar_zenith, &viewing_zenith, &scattering_angle, &albedo, &do_upwelling, temperature_layers.dataFirst(), air_number_density.dataFirst(), nm_grid_f.dataFirst(), nm_grid_out_f.dataFirst(), solar_irradiance_f.dataFirst(), total_optical_depth_f.dataFirst(), raman_spec.dataFirst(), &problems);

    if (problems) {
        throw Exception("raman_sioris fortran code encountered an internal issue.");
    }

    // Reverse back to order of calling grid
    if (reverse_grid_out) {
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
/// We currently only work with grounds that have a spurr_brdf_type
/// of LAMBERTIAN. We need to get the approximate albedo for the
/// surface, and just have the logic in place for that. We could
/// support other types if we work out the logic for getting an albedo
/// estimate (in "evaluate_albedo" function).
//-----------------------------------------------------------------------

RamanSiorisEffect::RamanSiorisEffect
(const SpectralDomain& Solar_and_odepth_spec_domain,
 double scale_factor,
 int channel_index,
 const DoubleWithUnit& solar_zenith, 
 const DoubleWithUnit& observation_zenith, 
 const DoubleWithUnit& relative_azimuth,
 const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
 const boost::shared_ptr<SolarModel>& solar_model,
 const boost::shared_ptr<StateMapping> mapping,
 bool do_upwelling)
: SpectrumEffectImpBase(scale_factor, mapping),
  solar_and_odepth_spec_domain_(Solar_and_odepth_spec_domain),
  channel_index_(channel_index),
  do_upwelling_(do_upwelling),
  atmosphere_(atmosphere),
  solar_model_(solar_model)
{
    // Convert angles to degrees since these should not change
    solar_zenith_ = solar_zenith.convert(units::deg).value;
    obs_zenith_ = observation_zenith.convert(units::deg).value;
    relative_azimuth_ = relative_azimuth.convert(units::deg).value; // stored for use in cloning
    scattering_angle_ = scattering_angle(observation_zenith, solar_zenith, relative_azimuth).convert(units::deg).value;

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
    firstIndex i1;
    secondIndex i2;
    // Convert dry air number density to the necessary units of molecules/cm^2
    Array<double, 1> dry_air_density = atmosphere_->absorber_ptr()->total_air_number_density_layer(channel_index_).convert(Unit("cm^-2")).value.value();

    // Compute total optical depth
    Array<double, 1> wn_grid = solar_and_odepth_spec_domain_.wavenumber();
    Array<double, 2> total_optical_depth(wn_grid.rows(),
					 temperature_layers_.rows());

    for(int wn_idx = 0; wn_idx < wn_grid.rows(); wn_idx++) {
        total_optical_depth(wn_idx, ra) =
	  atmosphere_->optical_depth_wrt_rt(wn_grid(wn_idx),
					    channel_index_).value();
    }
    
    // Absolute magnitude of the solar spectrum doesn't matter (we
    // calculate a ratio, so magnitude cancels out). No need to
    // convert to a certain unit;
    Array<double, 1> solar_spectrum = solar_model_->
      solar_spectrum(solar_and_odepth_spec_domain_).spectral_range().data();

    double wn_middle = (wn_grid(0) + wn_grid(wn_grid.rows() - 1)) / 2;
    double albedo = evaluate_albedo(wn_middle, channel_index_);

    // Compute raman spectrum
    Array<double, 1> raman_spec =
       compute_raman_sioris(solar_zenith_, obs_zenith_,
         scattering_angle_, albedo, do_upwelling_,
         temperature_layers_, dry_air_density, solar_and_odepth_spec_domain_,
         Spec.spectral_domain(), solar_spectrum, total_optical_depth);

    // Load scale factor
    AutoDerivative<double> scale_factor = coefficient()(0);

    // Apply raman scattering scaling.
    ArrayAd<double, 1> spec_rad = Spec.spectral_range().data_ad();

    // Chain rule for derivatives. We could just convert to and from
    // auto_derivative to do this automatically, but it isn't too hard
    // to just write out the explicit rule, which is faster to run.
    if(!spec_rad.is_constant() && !scale_factor.is_constant() &&
       spec_rad.number_variable() != scale_factor.number_variable())
	  throw Exception("Need to have the same number of variables for Spec jacobian and scale_factor"); 
    if(!spec_rad.is_constant() || !scale_factor.is_constant()) {
      spec_rad.resize_number_variable(std::max(spec_rad.number_variable(),
					       scale_factor.number_variable()));
      spec_rad.jacobian() = spec_rad.jacobian()(i1, i2) *
	(raman_spec(i1) * scale_factor.value() + 1.0);
      if(!scale_factor.is_constant())
	spec_rad.jacobian() += spec_rad.value()(i1) * raman_spec(i1) *
	  scale_factor.gradient()(i2);
    }
    spec_rad.value() *= (raman_spec * scale_factor.value() + 1.0);
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
/// Come up with an estimate of the surface albedo. This doesn't need
/// to be too accurate, the results don't depend strongly on the
/// albedo value. We currently look at the surface_parameter for the
/// middle wavenumber. We only work with LAMBERTIAN type.
//-----------------------------------------------------------------------

double RamanSiorisEffect::evaluate_albedo(double wn, int cindex) const
{
  // For ground observations, skip the surface contribution
  if(!do_upwelling_)
    return 0;
  if(atmosphere_->ground()->spurr_brdf_type() != SpurrBrdfType::LAMBERTIAN) {
    Exception e;
    e << "RamanSiorisEffect only works with a ground that has a SpurrBrdfType of LAMBERTIAN.\n"
      << "  ground: " << *atmosphere_->ground() << "\n";
    throw e;
  }
  return atmosphere_->ground()->surface_parameter(wn, cindex)(0).value();
}

//-----------------------------------------------------------------------
/// Create a new copy of this class using the same internal state 
/// as it currently exists
//-----------------------------------------------------------------------

boost::shared_ptr<SpectrumEffect> RamanSiorisEffect::clone() const
{

    return boost::make_shared<RamanSiorisEffect>
      (solar_and_odepth_spec_domain_, coefficient().value()(0),
       channel_index_,
       DoubleWithUnit(solar_zenith_, units::deg),
       DoubleWithUnit(obs_zenith_, units::deg),
       DoubleWithUnit(relative_azimuth_, units::deg),
       atmosphere_,
       solar_model_,
       mapping,
       do_upwelling_);
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
     << "  Do upwelling:           " << do_upwelling_ << "\n";
}
