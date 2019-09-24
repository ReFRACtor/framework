%include "fp_common.i"
%{
#include "pca_rt.h"
#include "sub_state_vector_array.h"
%}

%base_import(radiative_transfer_fixed_stokes_coefficient)

%import "pca_optical_properties.i"
%import "pca_binning.i"
%import "pca_eigensolver.i"
%import "atmosphere_standard.i"

%fp_shared_ptr(FullPhysics::PCARt);

namespace FullPhysics {

%feature("notabstract") PCARt;

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

    virtual int number_stream() const;
    virtual int number_stokes() const;

    virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
    virtual ArrayAd<double, 2> stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const;

    const boost::shared_ptr<AtmosphereStandard>& atmosphere() const;

    const boost::shared_ptr<PCAOpticalPropertiesAtmosphere> optical_properties() const;
    const boost::shared_ptr<PCABinning> binning() const;
    const boost::shared_ptr<PCAEigenSolver> solver(const int bin_index);
  
    virtual void print(std::ostream& Os, bool Short_form = false) const;
};
}
