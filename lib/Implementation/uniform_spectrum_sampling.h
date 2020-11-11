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

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
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
