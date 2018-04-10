#ifndef UNIFORM_SPECTRUM_SAMPLING_H
#define UNIFORM_SPECTRUM_SAMPLING_H

#include "spectrum_sampling.h"
#include "fp_exception.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This spectrum sampling creates a grid that begins the nearest
  ils_half_width before the low resolution grid beginning
  point and uniformly increases along some spacing until one
  ils_half_width past the end of the lowe resolution grid.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class UniformSpectrumSampling : public SpectrumSampling {
public:
    UniformSpectrumSampling(const ArrayWithUnit<double, 1>& Spec_spacing)
    	: spec_spacing(Spec_spacing), SpectrumSampling(Spec_spacing.rows()) { }
   
    virtual ~UniformSpectrumSampling() { }
  
    virtual SpectralDomain spectral_domain(int spec_index,
                                           const SpectralDomain& Lowres_grid, 
                                           const DoubleWithUnit& Ils_half_width) const;

    virtual void print(std::ostream& Os) const;
   
private:
    ArrayWithUnit<double, 1> spec_spacing;
};
}
#endif
