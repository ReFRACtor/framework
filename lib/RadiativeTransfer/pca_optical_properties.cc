#include "pca_optical_properties.h"

#include <boost/progress.hpp>

using namespace FullPhysics;
using namespace blitz;

PCAOpticalPropertiesAtmosphere::PCAOpticalPropertiesAtmosphere(const boost::shared_ptr<AtmosphereOco>& atm, const SpectralDomain& spec_domain, int channel_index, std::string primary_absorber, bool show_progress)
: atmosphere(atm), channel_index_(channel_index), show_progress_(show_progress)
{
    wavenumber_.reference(spec_domain.wavenumber());
   
    gas_optical_depth_.reference(Array<double, 2>(atmosphere->number_layer(), wavenumber_.rows(), blitz::ColumnMajorArray<2>()));
    total_optical_depth_.reference(Array<double, 2>(gas_optical_depth_.shape(), blitz::ColumnMajorArray<2>()));
    single_scattering_albedo_.reference(Array<double, 2>(gas_optical_depth_.shape(), blitz::ColumnMajorArray<2>()));
    primary_gas_dominates_.reference(Array<int, 1>(wavenumber_.rows(), blitz::ColumnMajorArray<1>()));

    int num_intermediate_vars = 2 + (atmosphere->aerosol_ptr() ? atmosphere->aerosol_ptr()->number_particle() : 0);
    intermediate_.reference(Array<double, 3>(atmosphere->number_layer(), num_intermediate_vars, wavenumber_.rows(), blitz::ColumnMajorArray<3>()));

    lambertian = boost::dynamic_pointer_cast<GroundLambertian>(atmosphere->ground());

    if (lambertian) {
        surface_albedo_.reference(Array<double, 1>(wavenumber_.rows(), blitz::ColumnMajorArray<1>()));
    }

    primary_abs_index_ = atmosphere->absorber_ptr()->gas_index(primary_absorber);

    if (primary_abs_index_ < 0) {
        Exception err;
        err << "Could not find index of primary absorber named: " << primary_absorber;
        throw err;
    }

    compute_properties();
}

void PCAOpticalPropertiesAtmosphere::compute_properties()
{
    Range all = Range::all();
    firstIndex i1; secondIndex i2;

    boost::shared_ptr<boost::progress_display> progress;
    if (show_progress_) {
        progress.reset(new boost::progress_display(wavenumber_.rows(), *Logger::stream()));
    }

    // Gather optical depth and single scattering albedo for binning routines
    for(int wn_idx = 0; wn_idx < wavenumber_.rows(); ++wn_idx) {
        Array<double, 2> gas_od(atmosphere->absorber_ptr()->optical_depth_each_layer(wavenumber_(wn_idx), channel_index_).value());
        
        // Check to see if primary absorber dominates column gas absorptions
        Array<double, 1> gas_col_tot(sum(gas_od, i1));

        if ((gas_col_tot(primary_abs_index_)/sum(gas_col_tot)) > 0.75) {
            primary_gas_dominates_(wn_idx) = 1;
        } else {
            primary_gas_dominates_(wn_idx) = 0;
        }

        gas_optical_depth_(all, wn_idx) = sum(gas_od, i2);
        total_optical_depth_(all, wn_idx) = atmosphere->optical_depth_wrt_iv(wavenumber_(wn_idx), channel_index_).value();
        single_scattering_albedo_(all, wn_idx) = atmosphere->single_scattering_albedo_wrt_iv(wavenumber_(wn_idx), channel_index_).value();

        intermediate_(all, all, wn_idx) = atmosphere->intermediate_variable(wavenumber_(wn_idx), channel_index_).value();
    
        if (lambertian) {
            surface_albedo_(wn_idx) = lambertian->albedo(DoubleWithUnit(wavenumber_(wn_idx), units::inv_cm), channel_index_).value();
        }

        if (progress) *progress += 1;
    }

}
