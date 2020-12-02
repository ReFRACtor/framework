#include "scalar_with_unit.h"
#include "spectral_domain.h"

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Variation of convert_wave that also handles the units of
/// sample_index. 
//-----------------------------------------------------------------------
template<class T>
ScalarWithUnit<T> ScalarWithUnit<T>::convert_wave
(const Unit& R, 
 const SpectralDomain& Pixel_grid) const
{
  if(units.is_commensurate(units::sample_index)) {
    int ind = (int) round(value) - 1;
    range_check(ind, 0, Pixel_grid.data().rows());
    ScalarWithUnit<T> d(Pixel_grid.data()(ind), Pixel_grid.units());
    return d.convert_wave(R);
  } else {
    return convert_wave(R);
  }
}
