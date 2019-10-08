#include "pca_binning.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *binlims, double *gasdat, int *ncnt, int *index, int *bin);

    void create_bin_uvvswir_v4(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);

    void create_bin_uvvswir_v5(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *ncnt, int *index, int *bin);
}

PCABinning::PCABinning(const boost::shared_ptr<PCAOpticalProperties>& optical_properties, const Method bin_method, const int num_bins)
: opt_props_(optical_properties), bin_method_(bin_method), num_bins_(num_bins)
{
    compute_bins();
}

void PCABinning::compute_bins()
{

    // Call selected binning routine
    int nlayer = opt_props_->gas_optical_depth().rows();
    int ndat = opt_props_->primary_gas_dominates().rows();
    
    // NCNT in the fortran binning code is allocated as 0:NBIN, so allocated as 1+ number of bins to avoid
    // invalid memory writes
    Array<int, 1> num_points(num_bins_+1, ColumnMajorArray<1>());

    // Binning routine packs all indexes one after another
    Array<int, 1> indexes_packed(ndat, ColumnMajorArray<1>());
    
    // Unused value but required to complete interface
    Array<int, 1> bins(ndat, ColumnMajorArray<1>());
    Array<double, 1> binlims(ndat, ColumnMajorArray<1>());

    // Ensure that inputs are all column major arrays by copying values from the optical properties
    Array<double, 2> taug(nlayer, ndat, ColumnMajorArray<2>());
    Array<double, 2> tau_tot(nlayer, ndat, ColumnMajorArray<2>());
    Array<double, 2> omega(nlayer, ndat, ColumnMajorArray<2>());

    taug = opt_props_->gas_optical_depth();
    tau_tot = opt_props_->total_optical_depth();
    omega = opt_props_->single_scattering_albedo();
    
    // Call fortran binning routine
    switch(bin_method_) {
        case UVVSWIR_V3:
            if(num_bins_ != 9) {
                throw Exception("UVVSWIR PCA binning V3 method must use 9 bins");
            }

            create_bin_uvvswir_v3(
                &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
                binlims.dataFirst(),
                taug.dataFirst(), 
                num_points.dataFirst(), indexes_packed.dataFirst(), bins.dataFirst());
            break;
        case UVVSWIR_V4:
            if(num_bins_ != 11) {
                throw Exception("UVVSWIR PCA binning V4 method must use 11 bins");
            }

            create_bin_uvvswir_v4(
                &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
                taug.dataFirst(), 
                tau_tot.dataFirst(), 
                omega.dataFirst(), 
                opt_props_->primary_gas_dominates().dataFirst(), 
                num_points.dataFirst(), indexes_packed.dataFirst(), bins.dataFirst());
            break;
        case UVVSWIR_V5:
            create_bin_uvvswir_v5(
                &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
                taug.dataFirst(), 
                tau_tot.dataFirst(), 
                omega.dataFirst(), 
                num_points.dataFirst(), indexes_packed.dataFirst(), bins.dataFirst());
            break;
        default:
            Exception err;
            err << "Unknown PCA binning method: " << bin_method_;
            throw err;
            break;
    }

    // Unpack indexes
    int pack_start = 0;
    for(int bidx = 0; bidx < num_bins_; bidx++) {
        int npoints = num_points(bidx);
        Array<int, 1> curr_indexes(npoints);
        // Remove 1 to make zero based indexes
        curr_indexes = indexes_packed(Range(pack_start, pack_start + npoints-1)) - 1;;
        pack_start += npoints;
        bin_indexes_.push_back(curr_indexes);
    }

    // Copy over bin number of points to a correctly sized object
    num_bin_points_.reference(Array<int, 1>(num_bins_, ColumnMajorArray<1>()));
    num_bin_points_ = num_points(Range(0, num_bins_-1));

}
