#include "planck.h"

using namespace blitz;
using namespace FullPhysics;

// Inform C about the Fortran function we will use
extern "C" {
    void planckfunction(const double* wavenumber, const double* temperature, double* bbfunc, int* smallv, bool* fail, const int* message_len, char* message);
    void planckfunction_plus(const double* wavenumber, const double* temperature, double* bbfunc, double* deriv_bbfunc, int* smallv, bool* fail, const int* message_len, char* message);
}

/// Computes the black body function value using the planck function for a given wavenumber (cm^-1) and temperature (K).
/// Result value will be in W/m^2/sr/cm^-1
///
/// Does not compute the derivative.
double FullPhysics::planck(double wn, double temperature)
{
    double bbfunc;
    int smallv;
    bool fail;
    const int message_len = 50;
    char message[message_len];

    planckfunction(&wn, &temperature, &bbfunc, &smallv, &fail, &message_len, message);

    if (fail) {
        Exception err;
        err << "plankfunction_plus failed with error message: " << message;
        throw err;
    }

    return bbfunc;
}

/// Computes the black body function value using the planck function for a given wavenumber (cm^-1) and temperature (K).
/// Result value will be in W/m^2/sr/cm^-1
///
/// Computes the derivative of the black body radiation with respect to temperature. Multiplies this by
/// the passed in gradient for use in the return value.
AutoDerivative<double> FullPhysics::planck(double wn, double temperature, Array<double, 1>& gradient)
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

    Array<double, 1> out_gradient(gradient.rows());
    out_gradient = gradient * deriv_bbfunc;

    AutoDerivative<double> bbody(bbfunc, out_gradient);

    return bbody;
}

