#include "pca_binning.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void create_bin_uvvswir_v3(int *E_nlayers, int *E_ndat, int *E_maxbins, int *ndat, int *nlay, int *nbin, double *gasdat, double *taudp, double *omega, int *absflag, int *ncnt, int *index, int *bin);

    void pca_eigensolver(int *Max_Eofs, int *maxpoints, int *maxlayers, int *maxlayers21, int *n_Eofs, int *npoints, int *nlayers, int *nlayers2, int *nlayers21, double *taudp, double *omega, double *Atmosmean, double *Eofs, double *PrinComps, bool *fail, int *message_len, char *message, int *trace_len, char *trace);
    void pca_eigensolver_alb(int *Max_Eofs, int *maxpoints, int *maxlayers, int *maxlayers21, int *n_Eofs, int *npoints, int *nlayers, int *nlayers2, int *nlayers21, double *taudp, double *omega, double *albedo, double *Atmosmean, double *Albmean, double *Eofs, double *PrinComps, bool *fail, int *message_len, char *message, int *trace_len, char *trace);
}

PCABinning::PCABinning(const boost::shared_ptr<PCAOpticalProperties>& optical_properties, int num_bins, int num_eofs)
: opt_props_(optical_properties), num_bins_(num_bins), num_eofs_(num_eofs)
{
    compute_bins();
    compute_eofs();
}

void PCABinning::compute_bins()
{

    // Call selected binning routine
    int nlayer = opt_props_->gas_optical_depth().rows();
    int ndat = opt_props_->primary_gas_dominates().rows();
    
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

void PCABinning::compute_eofs()
{
    int nlayers = opt_props_->gas_optical_depth().rows();
    int nlayers2 = nlayers*2;
    int nlayers21 = nlayers2 + 1;

    Array<double, 2> atmosmean(nlayers, 2, blitz::ColumnMajorArray<2>());
    Array<double, 2> albmean(nlayers, 2, blitz::ColumnMajorArray<2>());;
    Array<double, 2> eofs(num_eofs_, nlayers21, blitz::ColumnMajorArray<2>());;

    bool fail;
    int message_len = 100;
    int trace_len = 100;
    blitz::Array<char, 1> message(message_len);
    blitz::Array<char, 1> trace(trace_len);

    for(int bidx = 0; bidx < num_bins_; bidx++) {
        int npoints = num_bin_points_(bidx);

        Array<double, 2> princomps(num_eofs_, npoints, blitz::ColumnMajorArray<2>());;

        if (opt_props_->surface_albedo().rows() > 0) {
            // Use eigen solver capable of utilizing albedo if those values are available
            pca_eigensolver_alb(
                    &num_eofs_, &npoints, &nlayers, &nlayers21, &num_eofs_, &npoints, &nlayers, &nlayers2, &nlayers21, 
                    opt_props_->total_optical_depth().dataFirst(), 
                    opt_props_->single_scattering_albedo().dataFirst(), 
                    opt_props_->surface_albedo().dataFirst(), 
                    atmosmean.dataFirst(), albmean.dataFirst(), eofs.dataFirst(), princomps.dataFirst(),
                    &fail, &message_len, message.dataFirst(), &trace_len, trace.dataFirst());
        } else {
            // No surface albedo available
            pca_eigensolver(
                    &num_eofs_, &npoints, &nlayers, &nlayers21, &num_eofs_, &npoints, &nlayers, &nlayers2, &nlayers21, 
                    opt_props_->total_optical_depth().dataFirst(), 
                    opt_props_->surface_albedo().dataFirst(), 
                    atmosmean.dataFirst(), eofs.dataFirst(), princomps.dataFirst(),
                    &fail, &message_len, message.dataFirst(), &trace_len, trace.dataFirst());
        }
    }

}
