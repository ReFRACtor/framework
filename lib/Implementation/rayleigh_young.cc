#include "rayleigh_young.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void rayleigh_wrap(const double* wn,
                       const double* depolar_fact,
                       const double* a,
                       const double* b,
                       double* rayleigh_cross_section);
    /* depolar_fact depolarization factor for air is 0.02790 (young, 1980) */
    /* a and b are wavelength dependence coefficients for the refractive index
         (allen, 1964) (note: wavelengths must be in microns)
         For Earth, a = 2.871e-04, b = 5.67e-03 */
}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
RayleighYoung::RayleighYoung(const boost::shared_ptr<Pressure>& Pres,
                             const std::vector<boost::shared_ptr<Altitude> >& Alt,
                             const Constant& C)
    : RayleighImpBase(Pres, Alt, C),
      a(C.rayleigh_a().value), b(C.rayleigh_b().value),
      depolar_fact(C.rayleigh_depolarization_factor())
{
}

//-----------------------------------------------------------------------
/// Calculate the rayleigh cross section for the given
/// wavenumber/wavelength.
//-----------------------------------------------------------------------

DoubleWithUnit RayleighYoung::cross_section(const DoubleWithUnit& W) const
{
    double wn = W.convert_wave(units::inv_cm).value;

    double rayleigh_cross_section;
    rayleigh_wrap(&wn, &depolar_fact, &a, &b, &rayleigh_cross_section);

    return DoubleWithUnit(rayleigh_cross_section, Unit("m^2"));
}
