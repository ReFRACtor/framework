%include "fp_common.i"

%{
#include "planck.h"
%}

%import "auto_derivative.i"

namespace FullPhysics {

double planck(double wn, double temperature);
AutoDerivative<double> planck(double wn, double temperature, blitz::Array<double, 1>& gradient);

}
