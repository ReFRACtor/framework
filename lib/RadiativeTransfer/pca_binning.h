#ifndef PCA_BINNING_H
#define PCA_BINNING_H

#include "atmosphere_oco.h"
#include "ground_lambertian.h"
#include "spectral_domain.h"

namespace FullPhysics {

/****************************************************************//**
  Computes atmospheric optical properties need repeatedly in PCA
  calculations.
 *******************************************************************/

class PCAOpticalProperties {
public:

    PCAOpticalProperties(const boost::shared_ptr<AtmosphereOco>& atm, const SpectralDomain& spec_domain, int channel_index, std::string primary_absorber, bool show_progress=true);
    virtual ~PCAOpticalProperties() = default;

    blitz::Array<double, 1> wavenumber() const { return wavenumber_; }
    blitz::Array<double, 2> gas_optical_depth() const { return gas_optical_depth_; }
    blitz::Array<double, 2> total_optical_depth() const { return total_optical_depth_; }
    blitz::Array<double, 2> single_scattering_albedo() const { return single_scattering_albedo_; }
    blitz::Array<int, 1> primary_gas_dominates() const { return primary_gas_dominates_; }
    blitz::Array<double, 1> surface_albedo() const { return surface_albedo_; }

private:

    void compute_properties();

    boost::shared_ptr<AtmosphereOco> atmosphere;
    boost::shared_ptr<GroundLambertian> lambertian;

    bool show_progress_;

    int channel_index_;
    int primary_abs_index_;

    blitz::Array<double, 1> wavenumber_;
    blitz::Array<double, 2> gas_optical_depth_;
    blitz::Array<double, 2> total_optical_depth_;
    blitz::Array<double, 2> single_scattering_albedo_;
    blitz::Array<int, 1> primary_gas_dominates_;
    blitz::Array<double, 1> surface_albedo_;
};

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
    void compute_eofs();

    boost::shared_ptr<PCAOpticalProperties> opt_props_;
    int num_bins_;
    int num_eofs_;

    blitz::Array<int, 1> num_bin_points_;
    std::vector<blitz::Array<int, 1> > bin_indexes_;
 
};

}

#endif
