#ifndef PLANCK_H
#define PLANCK_H

#include "auto_derivative.h"

// Computes black body function plus derivative using LIDORT get_planckfunction_plus
namespace FullPhysics {

AutoDerivative<double> planck(double wn, double temperature);

}

#endif
