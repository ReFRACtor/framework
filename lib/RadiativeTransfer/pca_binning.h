#ifndef PCA_BINNING_H
#define PCA_BINNING_H

#include <boost/shared_ptr.hpp>

#include "generic_object.h"
#include "optical_properties_wrt_rt.h"

namespace FullPhysics {

/****************************************************************//**
  Compute PCA binned optical properties.
 *******************************************************************/

class PCABinning : public virtual GenericObject {
public:
    // Defined here so it doesn't end up in the global namespace for SWIG
    enum Method {
        UVVSWIR_V3 = 3,
        UVVSWIR_V4 = 4,
        UVVSWIR_V5 = 5,
    };

    PCABinning(const std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> >& optical_properties, const Method bin_method, const int num_bins, const int primary_absorber_index);
    virtual ~PCABinning() = default;

    /// Number of spectral points in each bin
    const blitz::Array<int, 1> num_bin_points() const { return num_bin_points_; }

    /// Index of spectral points into the full array of points for each bin
    const std::vector<blitz::Array<int, 1> > bin_indexes() const { return bin_indexes_; }

private:
    void compute_bins();

    std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> > opt_props_;

    Method bin_method_;
    int num_bins_;

    int primary_abs_index_;

    blitz::Array<int, 1> num_bin_points_;
    std::vector<blitz::Array<int, 1> > bin_indexes_;
 
};

}

#endif
