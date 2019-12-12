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

    Range r_gas(packed_idx, packed_idx + num_gas - 1);
    ArrayAd<double, 2> gas_od( packed_properties.value()(ra, r_gas), packed_properties.jacobian()(ra, r_gas, ra) );
    packed_idx += num_gas;

    Range r_aer_ext(packed_idx, packed_idx + num_aerosol - 1);
    ArrayAd<double, 2> aerosol_ext_od( packed_properties.value()(ra, r_aer_ext), packed_properties.jacobian()(ra, r_aer_ext, ra) );
    packed_idx += num_aerosol;

    Range r_aer_sca(packed_idx, packed_idx + num_aerosol - 1);
    ArrayAd<double, 2> aerosol_sca_od( packed_properties.value()(ra, r_aer_sca), packed_properties.jacobian()(ra, r_aer_sca, ra) );

    DoubleWithUnit spectral_point(wavenumber, units::inv_cm);
    boost::shared_ptr<AerosolPhaseFunctionHelper> aer_pf_helper(new AerosolPhaseFunctionComputeHelper(spectral_point, aerosol));

    initialize_with_jacobians(rayleigh_od, gas_od, aerosol_ext_od, aerosol_sca_od, aer_pf_helper);
}

ArrayAd<double, 2> OpticalPropertiesLsi::pack(const boost::shared_ptr<OpticalProperties>& source_properties)
{
    Range ra = Range::all();

    ArrayAd<double, 2> packed_v(
        source_properties->number_layers(), 
        1 + source_properties->number_gas_particles() + 2*source_properties->number_aerosol_particles(),
        source_properties->intermediate_jacobian().depth());

    int packed_idx = 0;

    ArrayAd<double, 2> gas_od(source_properties->gas_optical_depth_per_particle());
    ArrayAd<double, 1> ray_od(source_properties->rayleigh_optical_depth());
    ArrayAd<double, 2> aer_ext(source_properties->aerosol_extinction_optical_depth_per_particle());
    ArrayAd<double, 2> aer_sca(source_properties->aerosol_scattering_optical_depth_per_particle());

    // ******** TODO ************
    // Multiply jacobians by intermediate_jacobian so that they are in reference to the input jacobians
    // This way when we unpack averaged values they can be converted to wrt the RT parameters
    // Is this getting over convoluted?

    packed_v.value()(ra, packed_idx) = ray_od.value();
    packed_v.jacobian()(ra, packed_idx, ra) = ray_od.jacobian();
    packed_idx += 1;

    for(int gas_idx = 0; gas_idx < gas_od.cols(); gas_idx++) {
        packed_v.value()(ra, packed_idx) = gas_od.value()(ra, gas_idx);
        packed_v.jacobian()(ra, packed_idx, ra) = gas_od.jacobian()(ra, gas_idx, ra);
    }

    for(int aer_idx = 0; aer_idx < aer_ext.cols(); aer_idx++) {
        packed_v.value()(ra, packed_idx) = aer_ext.value()(ra, aer_idx); 
        packed_v.jacobian()(ra, packed_idx, ra) = aer_ext.jacobian()(ra, aer_idx, ra); 
        packed_idx += 1;
    }

    for(int aer_idx = 0; aer_idx < aer_ext.cols(); aer_idx++) {
        packed_v.value()(ra, packed_idx) = aer_sca.value()(ra, aer_idx); 
        packed_v.jacobian()(ra, packed_idx, ra) = aer_sca.jacobian()(ra, aer_idx, ra); 
        packed_idx += 1;
    }

    return packed_v;
}
