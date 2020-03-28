#include "float_with_unit.h"
#include "spectral_domain.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"

std::string float_with_unit_unit_get(const FloatWithUnit& V)
{
  return V.units.name();
}

void float_with_unit_unit_set(FloatWithUnit& V, std::string& Unit_name)
{
  V.units = Unit(Unit_name);
}

REGISTER_LUA_CLASS(FloatWithUnit)
.def(luabind::constructor<float, const std::string&>())
.def_readwrite("value", &FloatWithUnit::value)
.property("units", 
          &float_with_unit_unit_get,
          &float_with_unit_unit_set)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Variation of convert_wave that also handles the units of
/// sample_index. 
//-----------------------------------------------------------------------
    
FloatWithUnit FloatWithUnit::convert_wave
(const Unit& R, 
 const SpectralDomain& Pixel_grid) const
{
  if(units.is_commensurate(units::sample_index)) {
    int ind = (int) round(value) - 1;
    range_check(ind, 0, Pixel_grid.data().rows());
    FloatWithUnit d(Pixel_grid.data()(ind), Pixel_grid.units());
    return d.convert_wave(R);
  } else {
    return convert_wave(R);
  }
}
