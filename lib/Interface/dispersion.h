#ifndef DISPERSION_H
#define DISPERSION_H
#include "sample_grid.h"

/****************************************************************//**
The Dispersion type alias exists so that the older name for 
SampleGrid is still available to legacy code.
*******************************************************************/


namespace FullPhysics {
using Dispersion = SampleGrid;
}
#endif
