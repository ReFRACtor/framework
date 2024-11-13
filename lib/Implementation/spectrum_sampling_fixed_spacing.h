#ifndef SPECTRUM_SAMPLING_FIXED_SPACING_H
#define SPECTRUM_SAMPLING_FIXED_SPACING_H

#include "spectrum_sampling.h"

namespace FullPhysics {
/****************************************************************//**
  This generates a spectrum sampling that covers all the high
  resolution points needed to create the spectral domain covered by
  the given Instrument, subject to the given low resolution grid. For
  each range in the spectrum, we produce equally spaced points.

  This create high resolution points that are exact multiples of the
  spacing (so if the spacing is 0.01 you might have 2400.00 2400.01
  etc.).

  This is an ideal class to use if you have an Absco table (e.g.
  the forward model uses a AbscoAbsorber). By convention, the ABSCO
  tables are generated at exact multiples of the spacing. If we
  had some new Absco table where this wasn't the case, then this
  wouldn't be the spectrum sampling to use. We would need to create
  a new class that matched the Absco table (or alternatively allow
  interpolation of the Absco data to match this high resolution grid).

  This class is designed to *only* include high resolution grid points
  needed by the ILS with the given Edge_extension around the supplied
  Lowres_grid. So you don't get a uniform grid, in general there are
  high resolution points in a uniform grid that aren't include, For
  example the low resolution grid might be spaced more than 2 *
  Edge_extension, or we are missing low resolution grid point due to
  bad samples or something like that low This is generally what you
  want when using an ABSCO table.

  Contrast this with UniformSpectrumSampling, which really does
  include a uniform grid between a start and end point.
  
  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrumdoxygen for a description of each
  of these.
*******************************************************************/

class SpectrumSamplingFixedSpacing : public SpectrumSampling {
public:
  SpectrumSamplingFixedSpacing(const ArrayWithUnit<double, 1>& Spec_spacing)
    : SpectrumSampling(Spec_spacing.rows()), spec_spacing(Spec_spacing) { }

  virtual ~SpectrumSamplingFixedSpacing() {}
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Edge_extension) const;
  virtual void print(std::ostream& Os) const 
  { Os << "SpectrumSamplingFixedSpacing\n"
       << "  Spacing: " << spec_spacing << "\n";}
private:
  ArrayWithUnit<double, 1> spec_spacing;
  SpectrumSamplingFixedSpacing() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(SpectrumSamplingFixedSpacing);
#endif
