#include "pca_binning.h"

#include <boost/progress.hpp>

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);
}

void PCABinning::compute_bins(const SpectralDomain& spec_domain, int channel_index, std::string primary_absorber)
{
    Range all = Range::all();
    firstIndex i1; secondIndex i2;

    Array<double, 1> wavenumber(spec_domain.wavenumber());
    boost::shared_ptr<boost::progress_display> progress(new boost::progress_display(wavenumber.rows(), *Logger::stream()));
   
    Array<double, 2> gas_optical_depth(atmosphere->number_layer(), wavenumber.rows(), blitz::ColumnMajorArray<2>());
    Array<double, 2> total_optical_depth(gas_optical_depth.shape(), blitz::ColumnMajorArray<2>());
    Array<double, 2> single_scattering_albedo(gas_optical_depth.shape(), blitz::ColumnMajorArray<2>());
    Array<int, 1> primary_gas_dominates(wavenumber.rows(), blitz::ColumnMajorArray<1>()); 

    int primary_abs_index = atmosphere->absorber_ptr()->gas_index(primary_absorber);

    if (primary_abs_index < 0) {
        Exception err;
        err << "Could not find index of primary absorber named: " << primary_absorber;
        throw err;
    }

    // Gather optical depth and single scattering albedo for binning routines
    for(int wn_idx = 0; wn_idx < wavenumber.rows(); ++wn_idx) {
        Array<double, 2> gas_od(atmosphere->absorber_ptr()->optical_depth_each_layer(wavenumber(wn_idx), channel_index).value());
        
        // Check to see if primary absorber dominates column gas absorptions
        Array<double, 1> gas_col_tot(sum(gas_od, i1));

        if ((gas_col_tot(primary_abs_index)/sum(gas_col_tot)) > 0.75) {
            primary_gas_dominates(wn_idx) = 1;
        } else {
            primary_gas_dominates(wn_idx) = 0;
        }

        gas_optical_depth(all, wn_idx) = sum(gas_od, i2);
        total_optical_depth(all, wn_idx) = atmosphere->optical_depth_wrt_iv(wavenumber(wn_idx), channel_index).value();
        single_scattering_albedo(all, wn_idx) = atmosphere->single_scattering_albedo_wrt_iv(wavenumber(wn_idx), channel_index).value();

        if (progress) *progress += 1;
    }

    // Call selected binning routine
    int nlayer = atmosphere->number_layer();
    int ndat = wavenumber.rows();
    int nbin = 11;

    ncnt_.reference(Array<int, 1>(nbin+1, blitz::ColumnMajorArray<1>()));
    index_.reference(Array<int, 1>(wavenumber.rows(), blitz::ColumnMajorArray<1>()));
    bin_.reference(Array<int, 1>(wavenumber.rows(), blitz::ColumnMajorArray<1>()));
    
    create_bin_uvvswir_v3(&nlayer, &ndat, &nbin, &ndat, &nlayer, &nbin, gas_optical_depth.dataFirst(), total_optical_depth.dataFirst(), single_scattering_albedo.dataFirst(), primary_gas_dominates.dataFirst(), ncnt_.dataFirst(), index_.dataFirst(), bin_.dataFirst());

    // Compute binned optical properties, packing into the intermediate variable
}
