#ifndef SAMPLE_GRID_SPECTRAL_DOMAIN_H
#define SAMPLE_GRID_SPECTRAL_DOMAIN_H
#include "sample_grid.h"
#include "spectral_domain.h"
#include "unit.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This is an implementation of SampleGrid that uses a
  SpectralDomain to store the wavenumbers.

  Note that there are two minor variations of the sample grid
  spectral domain. The first wavenumber returned can either be the
  spectral domain at the value of "1", or a value of "0". By
  convention, the index is 1 based for GOSAT and OCO, but 0 based
  for FTS.
*******************************************************************/
class SampleGridSpectralDomain: public SubStateVectorArray<SampleGrid> {
public:
  SampleGridSpectralDomain(const SpectralDomain& Spec_domain,
		       const std::string& Band_name,
		       int Number_pixel,
		       bool Is_one_based);

  virtual ~SampleGridSpectralDomain() {}

  virtual boost::shared_ptr<SampleGrid> clone() const;

  virtual std::string sub_state_identifier() const { return "sample_grid/" + band_name_; }

  virtual SpectralDomain sample_grid() const;
  virtual void print(std::ostream& Os) const;
private:
  void initialize();
  SpectralDomain spec_domain;
  std::string band_name_;
  bool is_one_based;
  /// This is an array like 1,2,3 ... number_sample.
  blitz::Array<double, 1> index_array;
  // Very similar to index_array, but always 1 based.
  blitz::Array<int, 1> spectral_index;
};
}
#endif
