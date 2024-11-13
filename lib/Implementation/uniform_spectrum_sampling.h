#ifndef UNIFORM_SPECTRUM_SAMPLING_H
#define UNIFORM_SPECTRUM_SAMPLING_H

#include "spectrum_sampling.h"
#include "fp_exception.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This spectrum sampling creates a grid that begins the nearest
  edge extension amount before the low resolution grid beginning
  point and uniformly increases along some spacing until one
  unit of the edge extension past the end of the low resolution grid.

  This create high resolution points that are exact multiples of the
  spacing (so if the spacing is 0.01 you might have 2400.00 2400.01
  etc.).

  This class is strictly a uniform grid. We determine the lowest and
  highest high resolution grid point to cover the Lowres_grid, and the
  have high resolution grid points at all the multiples of the
  spacing. Note that in general, this may include more high resolution
  grid point than needed to create the low resolution grid using an
  ILS.

  Contrast this with SpectrumSamplingFixedSpacing, which only includes
  points needed by the ILS to create the low resolution
  grid. Generally if you are using Absco tables this is *not* the
  class you want, instead you want to use
  SpectrumSamplingFixedSpacing. But if you really do want a complete
  high resolution grid (e.g, you want to plot it, or there is some
  other reason you need all the points), then this is the class you
  want. 
  
  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrumdoxygen for a description of each
  of these.
*******************************************************************/
class UniformSpectrumSampling : public SpectrumSampling {
public:
  UniformSpectrumSampling(const ArrayWithUnit<double, 1>& Spec_spacing)
    : SpectrumSampling(Spec_spacing.rows()), spec_spacing(Spec_spacing) { }
   
  virtual ~UniformSpectrumSampling() { }
  
  virtual SpectralDomain spectral_domain(int spec_index,
					 const SpectralDomain& Lowres_grid, 
					 const DoubleWithUnit& Edge_extension) const;

  virtual void print(std::ostream& Os) const;
   
private:
  ArrayWithUnit<double, 1> spec_spacing;
  UniformSpectrumSampling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(UniformSpectrumSampling);
#endif
