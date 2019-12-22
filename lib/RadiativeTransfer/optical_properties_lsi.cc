#include "optical_properties_lsi.h"

using namespace FullPhysics;
using namespace blitz;

OpticalPropertiesLsi::OpticalPropertiesLsi(const ArrayAd<double, 2>& packed_properties, double wavenumber, const boost::shared_ptr<Aerosol>& aerosol, int num_gas, int num_aerosol)
: OpticalPropertiesImpBase()
{
    Range ra = Range::all();
    
    int packed_idx = 0;
    ArrayAd<double, 1> rayleigh_od( packed_properties.value()(ra, packed_idx), packed_properties.jacobian()(ra, packed_idx, ra) );
    packed_idx += 1;
    rayleigh_optical_depth_.reference(rayleigh_od);

    Range r_gas(packed_idx, packed_idx + num_gas - 1);
    ArrayAd<double, 2> gas_od( packed_properties.value()(ra, r_gas), packed_properties.jacobian()(ra, r_gas, ra) );
    packed_idx += num_gas;
    gas_optical_depth_per_particle_.reference(gas_od);

    Range r_aer_ext(packed_idx, packed_idx + num_aerosol - 1);
    ArrayAd<double, 2> aerosol_ext_od( packed_properties.value()(ra, r_aer_ext), packed_properties.jacobian()(ra, r_aer_ext, ra) );
    packed_idx += num_aerosol;
    aerosol_extinction_optical_depth_per_particle_.reference(aerosol_ext_od);

    Range r_aer_sca(packed_idx, packed_idx + num_aerosol - 1);
    ArrayAd<double, 2> aerosol_sca_od( packed_properties.value()(ra, r_aer_sca), packed_properties.jacobian()(ra, r_aer_sca, ra) );
    aerosol_scattering_optical_depth_per_particle_.reference(aerosol_sca_od);

    DoubleWithUnit spectral_point(wavenumber, units::inv_cm);
    aerosol_phase_function_helper_.reset(new AerosolPhaseFunctionComputeHelper(spectral_point, aerosol));

    intermediate_jacobian_.resize(rayleigh_optical_depth_.number_variable(), rayleigh_optical_depth_.number_variable());

    intermediate_jacobian_ = 0.0;
    for(int lay_idx = 0; lay_idx < rayleigh_optical_depth_.rows(); lay_idx++) {
        for(int jac_idx = 0; jac_idx < rayleigh_optical_depth_.number_variable(); jac_idx++) {
            intermediate_jacobian_(lay_idx, jac_idx, jac_idx) = 1.0;
        }
    }

    assert_sizes();

    cached_num_moments = -1;
    cached_num_scattering = -1;

    initialized = true;

}

ArrayAd<double, 2> OpticalPropertiesLsi::pack(const boost::shared_ptr<OpticalProperties>& source_properties)
{
    Range ra = Range::all();

    ArrayAd<double, 2> packed_v(
        source_properties->number_layers(), 
        1 + source_properties->number_gas_particles() + 2*source_properties->number_aerosol_particles(),
        source_properties->intermediate_jacobian().depth());
    packed_v = 0;

    ArrayAd<double, 1> ray_od(source_properties->rayleigh_optical_depth());
    ArrayAd<double, 2> gas_od(source_properties->gas_optical_depth_per_particle());
    ArrayAd<double, 2> aer_ext(source_properties->aerosol_extinction_optical_depth_per_particle());
    ArrayAd<double, 2> aer_sca(source_properties->aerosol_scattering_optical_depth_per_particle());

    int packed_idx = 0;

    packed_v.value()(ra, packed_idx) = ray_od.value();
    packed_v.jacobian()(ra, packed_idx, ra) = ray_od.jacobian();
    packed_idx += 1;

    Range r_gas_jac(0, gas_od.number_variable() - 1);
    for(int gas_idx = 0; gas_idx < gas_od.cols(); gas_idx++) {
        packed_v.value()(ra, packed_idx) = gas_od.value()(ra, gas_idx);
        packed_v.jacobian()(ra, packed_idx, r_gas_jac) = gas_od.jacobian()(ra, gas_idx, r_gas_jac);
        packed_idx += 1;
    }

    Range r_aer_ext_jac(0, aer_ext.number_variable() - 1);
    for(int aer_idx = 0; aer_idx < aer_ext.cols(); aer_idx++) {
        packed_v.value()(ra, packed_idx) = aer_ext.value()(ra, aer_idx); 
        packed_v.jacobian()(ra, packed_idx, r_aer_ext_jac) = aer_ext.jacobian()(ra, aer_idx, r_aer_ext_jac); 
        packed_idx += 1;
    }

    Range r_aer_sca_jac(0, aer_ext.number_variable() - 1);
    for(int aer_idx = 0; aer_idx < aer_ext.cols(); aer_idx++) {
        packed_v.value()(ra, packed_idx) = aer_sca.value()(ra, aer_idx); 
        packed_v.jacobian()(ra, packed_idx, r_aer_sca_jac) = aer_sca.jacobian()(ra, aer_idx, r_aer_sca_jac); 
        packed_idx += 1;
    }

    return packed_v;
}
