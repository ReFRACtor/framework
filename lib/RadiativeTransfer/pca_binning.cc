#include "pca_binning.h"

#include <boost/progress.hpp>

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);
}

PCAOpticalProperties::PCAOpticalProperties(const boost::shared_ptr<AtmosphereOco>& atm, const SpectralDomain& spec_domain, int channel_index, std::string primary_absorber, bool show_progress)
: atmosphere(atm), channel_index_(channel_index), show_progress_(show_progress)
{
    wavenumber_.reference(spec_domain.wavenumber());
   
    gas_optical_depth_.reference(Array<double, 2>(atmosphere->number_layer(), wavenumber_.rows(), blitz::ColumnMajorArray<2>()));
    total_optical_depth_.reference(Array<double, 2>(gas_optical_depth_.shape(), blitz::ColumnMajorArray<2>()));
    single_scattering_albedo_.reference(Array<double, 2>(gas_optical_depth_.shape(), blitz::ColumnMajorArray<2>()));
    primary_gas_dominates_.reference(Array<int, 1>(wavenumber_.rows(), blitz::ColumnMajorArray<1>()));

    primary_abs_index_ = atmosphere->absorber_ptr()->gas_index(primary_absorber);

    if (primary_abs_index_ < 0) {
        Exception err;
        err << "Could not find index of primary absorber named: " << primary_absorber;
        throw err;
    }

    compute_properties();
}

void PCAOpticalProperties::compute_properties()
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

        if (progress) *progress += 1;
    }

}

///////////////////////////////////////

PCABinning::PCABinning(const boost::shared_ptr<PCAOpticalProperties>& optical_properties, int num_bins)
: opt_props_(optical_properties), num_bins_(num_bins)
{
    compute_bins();
}

void PCABinning::compute_bins()
{

    // Call selected binning routine
    int nlayer = opt_props_->gas_optical_depth().rows();
    int ndat = opt_props_->wavenumber().rows();
    
    blitz::Array<int, 1> bin_;

    num_bin_points_.reference(Array<int, 1>(num_bins_, blitz::ColumnMajorArray<1>()));

    // Binning routine packs all indexes one after another
    Array<int, 1> indexes_packed(Array<int, 1>(ndat, blitz::ColumnMajorArray<1>()));
    
    // Unused value but required to complete interface
    Array<int, 1> bins(ndat, blitz::ColumnMajorArray<1>());
    
    // Call fortran binning routine
    create_bin_uvvswir_v3(
        &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
        opt_props_->gas_optical_depth().dataFirst(), 
        opt_props_->total_optical_depth().dataFirst(), 
        opt_props_->single_scattering_albedo().dataFirst(), 
        opt_props_->primary_gas_dominates().dataFirst(), 
        num_bin_points_.dataFirst(), indexes_packed.dataFirst(), bins.dataFirst());

    // Unpack indexes
    int pack_start = 0;
    for(int bidx = 0; bidx < num_bins_; bidx++) {
        int npoints = num_bin_points_(bidx);
        Array<int, 1> curr_indexes(npoints);
        curr_indexes = indexes_packed(Range(pack_start, pack_start + npoints));
        pack_start += npoints;
        bin_indexes_.push_back(curr_indexes);
    }

}
