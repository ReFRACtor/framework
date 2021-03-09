#ifndef SAMPLE_GRID_SPECTRAL_DOMAIN_H
#define SAMPLE_GRID_SPECTRAL_DOMAIN_H
#include "sample_grid_imp_base.h"
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
class SampleGridSpectralDomain: public SampleGridImpBase {
public:
    SampleGridSpectralDomain(const SpectralDomain& Spec_domain, const std::string& Band_name);

    virtual ~SampleGridSpectralDomain() = default;

    virtual boost::shared_ptr<SampleGrid> clone() const;

    virtual std::string sub_state_identifier() const
    {
        return "sample_grid/" + band_name_;
    }

    virtual SpectralDomain sample_grid() const;
    virtual void print(std::ostream& Os) const;
private:
    SpectralDomain spec_domain;
    std::string band_name_;
};
}
#endif
