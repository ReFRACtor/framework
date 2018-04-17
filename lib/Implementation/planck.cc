#include "planck.h"

using namespace FullPhysics;

// Inform C about the Fortran function we will use
extern "C" {
    void planckfunction_plus(const double* wavenumber, const double* temperature, double* bbfunc, double* deriv_bbfunc, int* smallv, bool* fail, const int* message_len, char* message);
}

/// Computes the black body function value using the planck function for a given wavenumber (cm^-1) and temperature (K).
/// Result value will be in W/m^2/sr/cm^-1
AutoDerivative<double> FullPhysics::planck(double wn, double temperature)
{
    double bbfunc;
    double deriv_bbfunc;
    int smallv;
    bool fail;
    const int message_len = 50;
    char message[message_len];

    planckfunction_plus(&wn, &temperature, &bbfunc, &deriv_bbfunc, &smallv, &fail, &message_len, message);

    if (fail) {
        Exception err;
        err << "plankfunction_plus failed with error message: " << message;
        throw err;
    }

    AutoDerivative<double> bbody(bbfunc, 0, 1);
    bbody.gradient() = deriv_bbfunc;

    return bbody;
}
