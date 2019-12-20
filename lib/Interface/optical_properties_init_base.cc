#include "optical_properties_init_base.h"

using namespace blitz;
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Initializes an OpticalProperties object using precomputed values:
/// In the parameters \f$l\f$ stands for layer index, \f$p\f$ stands for particle index.
/// \param rayleigh_od Rayleigh optical depth for each layer (\f$\tau_{ray,l}\f$)
/// \param gas_od Gas Absorber optical depth for each layer and gas particle type (\f$\tau_{gas,lp}\f$)
/// \param aerosol_ext_od Aerosol extinction optical depth for each layer and aerosol particle type (\f$\tau_{aer\_ext,lp}\f$)
/// \param aerosol_sca_od Aerosol scattering optical depth for each layer and aerosol particle type (\f$\tau_{aer\_sca,lp}\f$)
//-----------------------------------------------------------------------

void OpticalPropertiesInitBase::initialize(const ArrayAd<double, 1>& rayleigh_od, 
                                   const ArrayAd<double, 2>& gas_od,
                                   const ArrayAd<double, 2>& aerosol_ext_od,
                                   const ArrayAd<double, 2>& aerosol_sca_od,
                                   const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments)
{
    boost::shared_ptr<AerosolPhaseFunctionHelper> aer_pf_helper(new AerosolPhaseFunctionPassThruHelper(aerosol_pf_moments));

    initialize_with_jacobians(rayleigh_od, gas_od, aerosol_ext_od, aerosol_sca_od, aer_pf_helper);

    assert_sizes();

    cached_num_moments = -1;
    cached_num_scattering = -1;

    initialized = true;
}

//-----------------------------------------------------------------------
/// Initializes an OpticalProperties object from atmospheric objects
/// at a specific spectral point for a given instrument channel.
//-----------------------------------------------------------------------

void OpticalPropertiesInitBase::initialize(const DoubleWithUnit spectral_point,
                                   const int channel_index,
                                   const boost::shared_ptr<Absorber>& absorber,
                                   const boost::shared_ptr<Rayleigh>& rayleigh,
                                   const boost::shared_ptr<Aerosol>& aerosol)
{

    Range ra = Range::all();
    double wn = spectral_point.convert_wave(units::inv_cm).value;

    ArrayAd<double, 1> rayleigh_od(rayleigh->optical_depth_each_layer(wn, channel_index));
    ArrayAd<double, 2> gas_od(absorber->optical_depth_each_layer(wn, channel_index));

    ArrayAd<double, 2> aerosol_ext_od;
    ArrayAd<double, 2> aerosol_sca_od;

    boost::shared_ptr<AerosolPhaseFunctionHelper> aer_pf_helper;

    if (aerosol) {
        aerosol_ext_od.reference(aerosol->extinction_optical_depth_each_layer(wn));
        aerosol_sca_od.reference(aerosol->scattering_optical_depth_each_layer(wn));

        aer_pf_helper.reset(new AerosolPhaseFunctionComputeHelper(spectral_point, aerosol));
    }

    initialize_with_jacobians(rayleigh_od, gas_od, aerosol_ext_od, aerosol_sca_od, aer_pf_helper);

    assert_sizes();

    cached_num_moments = -1;
    cached_num_scattering = -1;

    initialized = true;
}
