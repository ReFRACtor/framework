#include "optical_properties_pca.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void OpticalPropertiesPca::serialize(Archive & ar,
				     const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(OpticalPropertiesImpBase);
}

FP_IMPLEMENT(OpticalPropertiesPca);
#endif

//-----------------------------------------------------------------------
/// Reconstruct a OpticalProperties object from packed optical properties
/// along with some additional information.
///
/// The packed properties do not include the per gas particle optical
/// depth values. An error will occur if it is queried. Only the total
/// gas optical depth is available.
//-----------------------------------------------------------------------

OpticalPropertiesPca::OpticalPropertiesPca
(const ArrayAd<double, 2>& packed_properties,
 double wavenumber,
 const boost::shared_ptr<AerosolOptical>& aerosol,
 int UNUSED(num_gas),
 int num_aerosol)
: OpticalPropertiesImpBase()
{
    Range ra = Range::all();

    intermediate_jacobian_.reference(packed_properties.jacobian());
    int num_rt_var = intermediate_jacobian_.cols();

    if(!packed_properties.is_constant() && num_rt_var != num_aerosol + 2)  {
        Exception err;
        err << "Packed properties jacobian size is incorrect with " << num_rt_var << " variables. " 
            << "Expected " << (num_aerosol + 2) << " variables.";
        throw err;
    }

    if(num_aerosol > 0 and !aerosol) {
        throw Exception("Number of aerosols is non zero but aerosol object is null");
    }

    int packed_idx = 0;
    gas_optical_depth_per_layer_.value().reference(packed_properties.value()(ra, packed_idx).copy());
    gas_optical_depth_per_layer_.resize_number_variable(num_rt_var);
    if(!packed_properties.is_constant()) {
        gas_optical_depth_per_layer_.jacobian()(ra, packed_idx) = 1.0;
    }
    packed_idx += 1;

    rayleigh_optical_depth_.value().reference(packed_properties.value()(ra, packed_idx).copy());
    rayleigh_optical_depth_.resize_number_variable(num_rt_var);
    if(!packed_properties.is_constant()) {
        rayleigh_optical_depth_.jacobian()(ra, packed_idx) = 1.0;
    }
    packed_idx += 1;

    if (aerosol) {
        Range r_aer_ext(packed_idx, packed_idx + num_aerosol - 1);
        aerosol_extinction_optical_depth_per_particle_.value().reference(packed_properties.value()(ra, r_aer_ext));
        aerosol_extinction_optical_depth_per_particle_.resize_number_variable(num_rt_var);
        if(!packed_properties.is_constant()) {
            for(int aer_idx = 0; aer_idx < num_aerosol; aer_idx++) {
                aerosol_extinction_optical_depth_per_particle_.jacobian()(ra, aer_idx, packed_idx) = 1.0;
                packed_idx += 1;
            }
        } else {
            packed_idx += num_aerosol;
        }

        // Compute aerosol scattering from aerosol optical depth using aerosol properties
        // like how this has historically been computed
        compute_aerosol_scattering(wavenumber, aerosol);
    }

    DoubleWithUnit spectral_point(wavenumber, units::inv_cm);
    aerosol_phase_function_helper_.reset(new AerosolPhaseFunctionComputeHelper(spectral_point, aerosol));

    assert_sizes();

    cached_num_moments = -1;
    cached_num_scattering = -1;

    initialized = true;

}

//-----------------------------------------------------------------------
/// Pack optical properties into a form that can be manipulated by PCA
//-----------------------------------------------------------------------

ArrayAd<double, 2> OpticalPropertiesPca::pack(const boost::shared_ptr<OpticalProperties>& source_properties)
{
    Range ra = Range::all();

    ArrayAd<double, 2> packed_v(
        source_properties->number_layers(), 
        2 + source_properties->number_aerosol_particles(),
        source_properties->intermediate_jacobian().depth());
    packed_v = 0;

    ArrayAd<double, 1> gas_od(source_properties->gas_optical_depth_per_layer());
    ArrayAd<double, 1> ray_od(source_properties->rayleigh_optical_depth());
    ArrayAd<double, 2> aer_ext(source_properties->aerosol_extinction_optical_depth_per_particle());

    if(!source_properties->is_constant())
      packed_v.jacobian() = source_properties->intermediate_jacobian();

    int packed_idx = 0;

    if(source_properties->number_gas_particles() > 0)
      packed_v.value()(ra, packed_idx) = gas_od.value();
    packed_idx += 1;

    packed_v.value()(ra, packed_idx) = ray_od.value();
    packed_idx += 1;

    if (source_properties->number_aerosol_particles() > 0) {
        Range r_aer(packed_idx, packed_idx + source_properties->number_aerosol_particles() - 1);
        packed_v.value()(ra, r_aer) = aer_ext.value();
    }

    return packed_v;
}

//-----------------------------------------------------------------------
/// Compute aerosol scattering value from aerosol extinction by way of
/// the aerosol single scattering albedo defined by the ratio of the
/// extinction and scattering aerosol properties.
//-----------------------------------------------------------------------

void OpticalPropertiesPca::compute_aerosol_scattering(double wavenumber, const boost::shared_ptr<AerosolOptical>& aerosol)
{
    firstIndex i1; secondIndex i2; thirdIndex i3;
    Range ra = Range::all();

    aerosol_scattering_optical_depth_per_particle_.reference(aerosol_extinction_optical_depth_per_particle_.copy());

    for(int particle_index = 0; particle_index < aerosol_extinction_optical_depth_per_particle_.cols(); particle_index++) {
        ArrayAd<double, 1> t = aerosol->aerosol_property(particle_index)->scattering_coefficient_each_layer(wavenumber);
        ArrayAd<double, 1> t2 = aerosol->aerosol_property(particle_index)->extinction_coefficient_each_layer(wavenumber);

        t.value() /= t2.value();

        if(!t.is_constant() && !aerosol_scattering_optical_depth_per_particle_.is_constant() && aerosol_scattering_optical_depth_per_particle_.number_variable() != t.number_variable()) {
            throw Exception("We don't currently have the code working correctly for combining intermediates and state vector derivatives. We'll need to think through how to do this.");
        }

        if(!t.is_constant()) {
            t.jacobian() = 1 / t2.value()(i1) * t.jacobian()(i1,i2) -
                t.value()(i1)/(t2.value()(i1) * t2.value()(i1))  * t2.jacobian()(i1,i2);
        }

        Array<double, 1> v(aerosol_scattering_optical_depth_per_particle_.value()(ra, particle_index));
        Array<double, 2> jac(aerosol_scattering_optical_depth_per_particle_.jacobian()(ra, particle_index, ra));
        v *= t.value();

        if(t.is_constant()) {
            jac = t.value()(i1) * jac(i1, i2);
        } else {
            jac = t.value()(i1) * jac(i1,i2) + v(i1) * t.jacobian()(i1,i2);
        }
    }
}

ArrayAd<double, 2> OpticalPropertiesPca::gas_optical_depth_per_particle() const
{
    throw Exception("OpticalPropertiesPca does not preserver per gas optical depth information.");
}
