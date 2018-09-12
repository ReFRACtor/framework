#ifndef PCA_BINNING_H
#define PCA_BINNING_H

#include "pca_optical_properties.h"

namespace FullPhysics {

/****************************************************************//**
  Compute PCA binned optical properties.
 *******************************************************************/

class PCABinning {
public:
    PCABinning(const boost::shared_ptr<PCAOpticalProperties>& optical_properties, int num_bins, int num_eofs);
    virtual ~PCABinning() = default;

    /// Number of spectral points in each bin
    const blitz::Array<int, 1> num_bin_points() const { return num_bin_points_; }

    /// Index of spectral points into the full array of points for each bin
    const std::vector<blitz::Array<int, 1> > bin_indexes() const { return bin_indexes_; }

private:
    void compute_bins();

    boost::shared_ptr<PCAOpticalProperties> opt_props_;
    int num_bins_;
    int num_eofs_;

    blitz::Array<int, 1> num_bin_points_;
    std::vector<blitz::Array<int, 1> > bin_indexes_;
 
};

}

#endif
