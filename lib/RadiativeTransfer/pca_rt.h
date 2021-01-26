#ifndef PCA_RT_H
#define PCA_RT_H

#include "radiative_transfer_fixed_stokes_coefficient.h"

#include "atmosphere_standard.h"
#include "aerosol_optical.h"

#include "lidort_rt.h"
#include "twostream_rt.h"
#include "first_order_rt.h"

#include "pca_binning.h"
#include "optical_properties_pca.h"
#include "pca_eigensolver.h"

namespace FullPhysics {

/****************************************************************//**
 Simple heler container class to encapsulate the optical properties
 computed from the PCA Solver
*******************************************************************/

struct PCABinOpticalProperties {
    boost::shared_ptr<OpticalPropertiesPca> mean;
    std::vector<boost::shared_ptr<OpticalPropertiesPca> > eof_plus;
    std::vector<boost::shared_ptr<OpticalPropertiesPca> > eof_minus;
};

/****************************************************************//**
 Implements a radiative transfer using the PCA method where
 the requested spectral domain is broken up into bins based on
 their optical properties.
 
 Each bin has higher accuracy computations done for the mean value 
 and a certain number of EOF perturbations of the mean. A 
 correction  factor is determined how to combine the binned 
 computations with low accuracy RT computation.

 This class uses LIDORT, 2steam and First Order together to 
 speed up multiple scattering RT computations.
*******************************************************************/

class PCARt : public RadiativeTransferFixedStokesCoefficient {
public:
    PCARt(const boost::shared_ptr<AtmosphereStandard>& Atm,
          const std::string Primary_absorber,
          const PCABinning::Method Bin_method, const int Num_bins,
          const int Num_eofs,
          const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
          const blitz::Array<double, 1>& Sza, 
          const blitz::Array<double, 1>& Zen, 
          const blitz::Array<double, 1>& Azm,
          int Number_streams, 
          int Number_moments, 
          bool do_solar_sources = true, 
          bool do_thermal_emission = false);

    virtual ~PCARt() = default;

    virtual int number_stream() const { return lidort_rt->number_stream(); }

    virtual int number_stokes() const { return stokes_coef->stokes_coefficient().cols(); }

    virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
    virtual ArrayAd<double, 2> stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const;

    const boost::shared_ptr<AtmosphereStandard>& atmosphere() const { return atm; }

    const boost::shared_ptr<LidortRt> lidort() const { return lidort_rt; }
    const boost::shared_ptr<TwostreamRt> twostream() const { return twostream_rt; }
    const boost::shared_ptr<FirstOrderRt> first_order() const { return first_order_rt; }

    // Debugging accessors
    const std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> > optical_properties() const { return pca_opt; }
    const boost::shared_ptr<PCABinning> binning() const { return pca_bin; }
    const boost::shared_ptr<PCAEigenSolver> solver(const int bin_index) { 

        if (bin_index < 0 || bin_index >= pca_solvers.size()) {
            Exception err;
            err << "Index for binned eigen solver: " << bin_index << " exceeds size of solvers saved: " << pca_solvers.size();
            throw err;
        }
        return pca_solvers[bin_index];

    }

    virtual void print(std::ostream& Os, bool Short_form = false) const;

private:

    // Steps of the PCA process
    void clear_pca_objects() const;
    void compute_bins(const SpectralDomain& Spec_domain, int Spec_index) const;
    const boost::shared_ptr<PCAEigenSolver> compute_bin_solution(const blitz::Array<int, 1>& data_indexes) const;
    const boost::shared_ptr<PCABinOpticalProperties> compute_bin_optical_props(boost::shared_ptr<PCAEigenSolver> pca_solver, double bin_wn) const;
    const double bin_effective_wavenumber(blitz::Array<double, 1> &win_wavenumbers, int bin_index) const;
    blitz::Array<double, 2> compute_bin_correction_factors(boost::shared_ptr<PCAEigenSolver>& pca_solver, boost::shared_ptr<PCABinOpticalProperties>& bin_opt_props, double bin_wn, int channel_index) const;

    boost::shared_ptr<AtmosphereStandard> atm;

    std::string primary_absorber;

    PCABinning::Method bin_method;
    int num_bins;
    int num_eofs;

    // Sizes needed for allocating arrays
    int num_gas;
    int num_aerosol;
    int num_layer;
    int num_packed_var;

    // This more specific interface is needed vs the abstract Aerosol interface
    boost::shared_ptr<AerosolOptical> aerosol_optical;

    boost::shared_ptr<LidortRt> lidort_rt;
    boost::shared_ptr<TwostreamRt> twostream_rt;
    boost::shared_ptr<FirstOrderRt> first_order_rt;

    // These are stored for the current stokes call for debugging purposes
    // They are empty until stokes or stokes_and_jacobian are called.
    // They are reset for each call
    mutable std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> > pca_opt;
    mutable boost::shared_ptr<PCABinning> pca_bin;
    mutable std::vector<boost::shared_ptr<PCAEigenSolver> > pca_solvers;
    
};
}
#endif
