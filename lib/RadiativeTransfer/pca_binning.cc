#include "pca_binning.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *binlims, double *gasdat, int *ncnt, int *index, int *bin);

    void create_bin_uvvswir_v4(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);

    void create_bin_uvvswir_v5(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *ncnt, int *index, int *bin);
}

PCABinning::PCABinning(const std::vector<boost::shared_ptr<OpticalProperties> >& optical_properties, const Method bin_method, const int num_bins, const int primary_absorber_index)
: opt_props_(optical_properties), bin_method_(bin_method), num_bins_(num_bins), primary_abs_index_(primary_absorber_index)
{
    compute_bins();

    // Check for consistency
    for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {
        if(bin_indexes_[bin_idx].rows() != num_bin_points_(bin_idx)) {
            Exception err;
            err << "Number of bin data indexes: " << bin_indexes_[bin_idx].rows()
                << " does not match number of bin points: " << num_bin_points_(bin_idx) 
                << " at bin " << bin_idx;
            throw err;
        }
    }

}

void PCABinning::compute_bins()
{
    int ndat = opt_props_.size();

    if (ndat == 0) {
        throw Exception("No optical property information present for PCA binning");
    }

    Range ra = Range::all();

    // Call selected binning routine
    int nlayer = opt_props_[0]->number_layers();
    
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

    for(int dom_idx = 0; dom_idx < ndat; dom_idx++) {
        taug(ra, dom_idx) = opt_props_[dom_idx]->gas_optical_depth_per_layer().value();
        tau_tot(ra, dom_idx) = opt_props_[dom_idx]->total_optical_depth().value();
        omega(ra, dom_idx) = opt_props_[dom_idx]->total_single_scattering_albedo().value();
    }

    // Only used for UVVSWIR_V4, resized there
    Array<int, 1> primary_gas_dominates = Array<int, 1>(blitz::ColumnMajorArray<1>());
    
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

            primary_gas_dominates.resize(ndat);

            for(int dom_idx = 0; dom_idx < ndat; dom_idx++) {
                Array<double, 2> gas_od(opt_props_[dom_idx]->gas_optical_depth_per_particle().value());
                Array<double, 1> gas_col_tot(gas_od.cols());
                for(int gas_idx = 0; gas_idx < gas_od.cols(); gas_idx++) {
                    gas_col_tot(gas_idx) = sum(gas_od(ra, gas_idx));
                }

                if ((gas_col_tot(primary_abs_index_)/sum(gas_col_tot)) > 0.75) {
                    primary_gas_dominates(dom_idx) = 1;
                } else {
                    primary_gas_dominates(dom_idx) = 0;
                }
            }

            create_bin_uvvswir_v4(
                &nlayer, &ndat, &num_bins_, &ndat, &nlayer, &num_bins_, 
                taug.dataFirst(), 
                tau_tot.dataFirst(), 
                omega.dataFirst(), 
                primary_gas_dominates.dataFirst(), 
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
        if (npoints > 0) {
            // Remove 1 to make zero based indexes
            curr_indexes = indexes_packed(Range(pack_start, pack_start + npoints-1)) - 1;;
            pack_start += npoints;
        }
        bin_indexes_.push_back(curr_indexes);
    }

    // Copy over bin number of points to a correctly sized object
    num_bin_points_.reference(Array<int, 1>(num_bins_, ColumnMajorArray<1>()));
    num_bin_points_ = num_points(Range(0, num_bins_-1));

}
