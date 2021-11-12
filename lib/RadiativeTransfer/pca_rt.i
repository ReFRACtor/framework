// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pca_rt.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(radiative_transfer_fixed_stokes_coefficient)

%import "atmosphere_standard.i"

%import "lidort_rt.i"
%import "twostream_rt.i"
%import "first_order_rt.i"

%import "optical_properties_wrt_rt.i"
%import "pca_binning.i"
%import "pca_eigensolver.i"

%fp_shared_ptr(FullPhysics::PCARt);

namespace FullPhysics {

%feature("notabstract") PCARt;

class PCARt : public RadiativeTransferFixedStokesCoefficient,
              public Observer<RtAtmosphere> {
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
	bool do_thermal_emission = false,
	bool do_3M_correction = false);

  %python_attribute(number_stream, virtual int);
  %python_attribute(number_stokes, virtual int);

  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain,
					 int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const SpectralDomain& Spec_domain, int Spec_index) const;

  %python_attribute(atmosphere,boost::shared_ptr<AtmosphereStandard>);
  %python_attribute(lidort, boost::shared_ptr<LidortRt>);
  %python_attribute(twostream, boost::shared_ptr<TwostreamRt>);
  %python_attribute(first_order, boost::shared_ptr<FirstOrderRt>);
  %python_attribute(optical_properties, std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> >);
  %python_attribute(binning, boost::shared_ptr<PCABinning>);
  const boost::shared_ptr<PCAEigenSolver> solver(const int bin_index);
  
  %pickle_serialization();
};
}
