#include "pca_binning.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);
}

PCABinning::PCABinning(const boost::shared_ptr<PCAOpticalProperties>& optical_properties, int num_bins)
: opt_props_(optical_properties), num_bins_(num_bins)
{
    compute_bins();
}

void PCABinning::compute_bins()
{

    // Call selected binning routine
    int nlayer = opt_props_->gas_optical_depth().rows();
    int ndat = opt_props_->primary_gas_dominates().rows();
    
    Array<int, 1> bin_;

    num_bin_points_.reference(Array<int, 1>(num_bins_, ColumnMajorArray<1>()));

    // Binning routine packs all indexes one after another
    Array<int, 1> indexes_packed(Array<int, 1>(ndat, ColumnMajorArray<1>()));
    
    // Unused value but required to complete interface
    Array<int, 1> bins(ndat, ColumnMajorArray<1>());

    // Ensure that inputs are all column major arrays by copying values from the optical properties
    Array<double, 2> taug(nlayer, ndat, ColumnMajorArray<2>());
    Array<double, 2> tau_tot(nlayer, ndat, ColumnMajorArray<2>());
    Array<double, 2> omega(nlayer, ndat, ColumnMajorArray<2>());

    taug = opt_props_->gas_optical_depth();
    tau_tot = opt_props_->total_optical_depth();
    omega = opt_props_->single_scattering_albedo();
    
    // Call fortran binning routine
    create_bin_uvvswir_v3(
        &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
        taug.dataFirst(), 
        tau_tot.dataFirst(), 
        omega.dataFirst(), 
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
