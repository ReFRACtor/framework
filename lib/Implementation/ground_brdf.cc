#include "ground_brdf.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

// Number of coefficients per band in the state vector
#define NUM_COEFF 7

// Number of parameters in spars for the RT
#define NUM_PARAMS 5

extern "C" {
    double black_sky_albedo_veg_f(const double* params, const double* sza);
    double black_sky_albedo_soil_f(const double* params, const double* sza);
    double exact_brdf_value_veg_f(const double* params, const double* sza, const double* vza, const double* azm);
    double exact_brdf_value_soil_f(const double* params, const double* sza, const double* vza, const double* azm);
}

#ifdef HAVE_LUA
#include "register_lua.h"

double black_sky_albedo_simple_veg(const blitz::Array<double, 1>& params, double sza) {
    return black_sky_albedo_veg_f(params.dataFirst(), &sza);
}

double black_sky_albedo_simple_soil(const blitz::Array<double, 1>& params, double sza) {
    return black_sky_albedo_soil_f(params.dataFirst(), &sza);
}

double exact_brdf_value_simple_veg(const blitz::Array<double, 1>& params, double sza, double vza, double azm) {
    return exact_brdf_value_veg_f(params.dataFirst(), &sza, &vza, &azm);
}

double exact_brdf_value_simple_soil(const blitz::Array<double, 1>& params, double sza, double vza, double azm) {
    return exact_brdf_value_soil_f(params.dataFirst(), &sza, &vza, &azm);
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfVeg, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_veg)
]
.scope
[
    luabind::def("kernel_value", &exact_brdf_value_simple_veg)
]
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBrdfSoil, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_soil)
]
.scope
[
    luabind::def("kernel_value", &exact_brdf_value_simple_soil)
]
REGISTER_LUA_END()
#endif

/****************************************************************//**
  Constructor that defines coefficients in a 2d array:
  Num_spectrometer * NUM_COEFF
  Each row has the NUM_COEFF BRDF parameters for a spectrometer.
  Coefficients are ordered:
  0: BRDF overall weight intercept
  1: BRDF overall weight slope
  2: Rahman kernel factor
  3: Rahman hotspot parameter
  4: Rahman asymmetry factor
  5: Rahman anisotropy parameter
  6: Breon kernel factor
 *******************************************************************/

GroundBrdf::GroundBrdf(const blitz::Array<double, 2>& Coeffs,
                       const blitz::Array<bool, 2>& Flag,
                       const ArrayWithUnit<double, 1>& Ref_points,
                       const std::vector<std::string>& Desc_band_names) 
: reference_points(Ref_points), desc_band_names(Desc_band_names)
{
    if(Coeffs.cols() != NUM_COEFF) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " is not " << NUM_COEFF << " as expected";
        throw err_msg;
    }

    if(Coeffs.rows() != Flag.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Flag: " << Flag.rows();
        throw err_msg;
    }

    if(Coeffs.cols() != Flag.cols()) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " does not match the number in Flag: " << Flag.cols();
        throw err_msg;
    }

    if(Coeffs.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    if(Ref_points.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Ref_points: " << Ref_points.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    // Make local arrays to deal with const issues on call to init. The init routine copies the data
    Array<double, 2> coeffs(Coeffs);
    Array<bool, 2> flags(Flag);

    // Flatten arrays for state vector
    init(Array<double, 1>(coeffs.dataFirst(), TinyVector<int, 1>(Coeffs.rows() * Coeffs.cols()), blitz::neverDeleteData),
         Array<bool, 1>(flags.dataFirst(), TinyVector<int, 1>(Flag.rows() * Flag.cols()), blitz::neverDeleteData));
}

/// Protected constructor that matches the dimensionality of coeff and flag arrays
GroundBrdf::GroundBrdf(const blitz::Array<double, 1>& Spec_coeffs,
                       const blitz::Array<bool, 1>& Flag, 
                       const ArrayWithUnit<double, 1>& Ref_points,
                       const std::vector<std::string>& Desc_band_names)
  : SubStateVectorArray<Ground>(Spec_coeffs, Flag),
    reference_points(Ref_points), desc_band_names(Desc_band_names)
{
}

ArrayAd<double, 1> GroundBrdf::surface_parameter(double wn, int spec_index) const
{
    AutoDerivative<double> w = weight(wn, spec_index);
    ArrayAd<double, 1> spars;
    spars.resize(NUM_PARAMS, coefficient().number_variable());
    spars(0) = w * rahman_factor(spec_index);
    spars(1) = hotspot_parameter(spec_index);
    spars(2) = asymmetry_parameter(spec_index);
    spars(3) = anisotropy_parameter(spec_index);
    spars(4) = w * breon_factor(spec_index);
    return spars;
}

const AutoDerivative<double> GroundBrdf::weight(double wn, int spec_index) const
{
    double ref_wn = reference_points(spec_index).convert_wave(units::inv_cm).value;
    ArrayAd<double, 1> weight_params(2, weight_intercept(spec_index).number_variable());
    weight_params(0) = weight_slope(spec_index);
    weight_params(1) = weight_intercept(spec_index);
    Poly1d weight_poly(weight_params, true);
    AutoDerivative<double> wn_ad(wn); // Make sure we use the AutoDerivative interface to Poly1d
    return weight_poly(wn_ad - ref_wn);
}

const AutoDerivative<double> GroundBrdf::weight_intercept(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + BRDF_WEIGHT_INTERCEPT_INDEX);
}

const AutoDerivative<double> GroundBrdf::weight_slope(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + BRDF_WEIGHT_SLOPE_INDEX);
}

const AutoDerivative<double> GroundBrdf::rahman_factor(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + RAHMAN_KERNEL_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::hotspot_parameter(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + RAHMAN_OVERALL_AMPLITUDE_INDEX);
}

const AutoDerivative<double> GroundBrdf::asymmetry_parameter(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + RAHMAN_ASYMMETRY_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::anisotropy_parameter(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + RAHMAN_GEOMETRIC_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::breon_factor(int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + BREON_KERNEL_FACTOR_INDEX);
}

//----

void GroundBrdf::weight_intercept(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + BRDF_WEIGHT_INTERCEPT_INDEX) = val;
}

void GroundBrdf::weight_slope(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + BRDF_WEIGHT_SLOPE_INDEX) = val;
}

void GroundBrdf::rahman_factor(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + RAHMAN_KERNEL_FACTOR_INDEX) = val;
}

void GroundBrdf::hotspot_parameter(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + RAHMAN_OVERALL_AMPLITUDE_INDEX) = val;
}

void GroundBrdf::asymmetry_parameter(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + RAHMAN_ASYMMETRY_FACTOR_INDEX) = val;
}

void GroundBrdf::anisotropy_parameter(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + RAHMAN_GEOMETRIC_FACTOR_INDEX) = val;
}

void GroundBrdf::breon_factor(int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + BREON_KERNEL_FACTOR_INDEX) = val;
}

const blitz::Array<double, 2> GroundBrdf::brdf_covariance(int spec_index) const
{ 
    range_check(spec_index, 0, number_spectrometer());

    // Return empty array if covariance is empty due to not being retrieved
    if(product(statevector_covariance().shape()) == 0) {
        return Array<double, 2>(0, 0);
    }

    int ret_num_param = statevector_covariance().rows() / desc_band_names.size();
    int ret_cov_offset = ret_num_param * spec_index;

    blitz::Array<double, 2> cov(NUM_COEFF, NUM_COEFF);
    cov = 0.0;

    // Copy out the retrieved covariance values into a matrix for the spectrometer that includes
    // zeros for non retrieved elements
    int in_idx_a = ret_cov_offset;
    for (int out_idx_a = 0; out_idx_a < NUM_COEFF; out_idx_a++) {
        if (used_flag_value()(out_idx_a)) {
            int in_idx_b = ret_cov_offset;
            for (int out_idx_b = 0; out_idx_b < NUM_COEFF; out_idx_b++) {
                if (used_flag_value()(out_idx_b)) {
                    cov(out_idx_a, out_idx_b) = statevector_covariance()(in_idx_a, in_idx_b++);
                }
            }
            in_idx_a++;
        }
    }

    return cov;
}

// Helper function 
blitz::Array<double, 1> GroundBrdf::black_sky_params(int Spec_index)
{
    double ref_wn = reference_points(Spec_index).convert_wave(units::inv_cm).value;
    double w = weight(ref_wn, Spec_index).value();

    blitz::Array<double, 1> params(NUM_PARAMS, blitz::ColumnMajorArray<1>());
    params(0) = w * rahman_factor(Spec_index).value();
    params(1) = hotspot_parameter(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = anisotropy_parameter(Spec_index).value();
    params(4) = w * breon_factor(Spec_index).value();
    return params;
}

// Helper function 
blitz::Array<double, 1> GroundBrdf::kernel_value_params(int Spec_index)
{
    blitz::Array<double, 1> params(NUM_PARAMS, blitz::ColumnMajorArray<1>());
    params(0) = rahman_factor(Spec_index).value();
    params(1) = hotspot_parameter(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = anisotropy_parameter(Spec_index).value();
    params(4) = breon_factor(Spec_index).value();
    return params;
}

double GroundBrdfVeg::black_sky_albedo(int Spec_index, double Sza)
{
    blitz::Array<double, 1> params = black_sky_params(Spec_index);
    return black_sky_albedo_veg_f(params.dataFirst(), &Sza);
}

double GroundBrdfSoil::black_sky_albedo(int Spec_index, double Sza)
{
    blitz::Array<double, 1> params = black_sky_params(Spec_index);
    return black_sky_albedo_soil_f(params.dataFirst(), &Sza);
}

double GroundBrdfVeg::kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm)
{
    return exact_brdf_value_veg_f(params.dataFirst(), &Sza, &Vza, &Azm);
}

double GroundBrdfVeg::kernel_value(int Spec_index, double Sza, double Vza, double Azm)
{
    blitz::Array<double, 1> params = kernel_value_params(Spec_index);
    return kernel_value_at_params(params, Sza, Vza, Azm);
}

double GroundBrdfSoil::kernel_value_at_params(const blitz::Array<double, 1>& params, double Sza, double Vza, double Azm)
{
    return exact_brdf_value_soil_f(params.dataFirst(), &Sza, &Vza, &Azm);
}

double GroundBrdfSoil::kernel_value(int Spec_index, double Sza, double Vza, const double Azm)
{
    blitz::Array<double, 1> params = kernel_value_params(Spec_index);
    return kernel_value_at_params(params, Sza, Vza, Azm);
}

std::string GroundBrdf::state_vector_name_i(int i) const {
    int b_idx = int(i / NUM_COEFF);
    int c_idx = i - NUM_COEFF * b_idx;

    std::stringstream name;
    name << "Ground BRDF " << breon_type() << " " << desc_band_names[b_idx] << " ";
    switch (c_idx) {
    case BRDF_WEIGHT_INTERCEPT_INDEX:
        name << "BRDF Weight Intercept";
        break;
    case BRDF_WEIGHT_SLOPE_INDEX:
        name << "BRDF Weight Slope";
        break;
    case RAHMAN_KERNEL_FACTOR_INDEX:
        name << "Rahman Factor";
        break;
    case RAHMAN_OVERALL_AMPLITUDE_INDEX:
        name << "Hotspot Parameter";
        break;
    case RAHMAN_ASYMMETRY_FACTOR_INDEX:
        name << "Asymmetry Parameter";
        break;
    case RAHMAN_GEOMETRIC_FACTOR_INDEX:
        name << "Anisotropy Parameter";
        break;
    case BREON_KERNEL_FACTOR_INDEX:
        name << "Breon Factor";
        break;
    default:
        name << "Unknown Index " << i;
    }

    return name.str();
}

void GroundBrdf::print(std::ostream& Os) const
{
    Os << "GroundBrdf:\n";
    for(int spec_idx = 0; spec_idx < number_spectrometer(); spec_idx++) {
        Os << "    " << desc_band_names[spec_idx] << ":" << std::endl;
        OstreamPad opad(Os, "        ");
        opad << "BRDF Weight Intercept: " << weight_intercept(spec_idx).value() << std::endl
             << "BRDF Weight Slope: " << weight_slope(spec_idx).value() << std::endl
             << "Rahman Factor: " << rahman_factor(spec_idx).value() << std::endl
             << "Hotspot Parameter: " << hotspot_parameter(spec_idx).value() << std::endl
             << "Asymmetry Parameter: " << asymmetry_parameter(spec_idx).value() << std::endl
             << "Anisotropy Parameter: " << anisotropy_parameter(spec_idx).value() << std::endl
             << "Breon Factor: " << breon_factor(spec_idx).value() << std::endl;
        opad.strict_sync();
    }
}
