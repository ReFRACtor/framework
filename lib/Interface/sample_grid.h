#ifndef SAMPLE_GRID_H
#define SAMPLE_GRID_H

#include "spectral_domain.h"

namespace FullPhysics {
/****************************************************************//**
  This determines the sample grid that should be used for each
  of the instrument indexes.

*******************************************************************/
class SampleGrid : public virtual GenericObject {
public:
  virtual ~SampleGrid() {}

//-----------------------------------------------------------------------
/// Wavenumbers/Wavelengths to use for the given spectrometer.
//-----------------------------------------------------------------------

  virtual SpectralDomain sample_grid(int Spec_index) const = 0;

protected:
//-----------------------------------------------------------------------
/// Default constructor
//-----------------------------------------------------------------------
  SampleGrid() {}
};
}
#endif
