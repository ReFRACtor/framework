#include "rayleigh_young.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(RayleighYoung, Rayleigh)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
                          const std::vector<boost::shared_ptr<Altitude> >&,
                          const boost::shared_ptr<Constant>&>())
REGISTER_LUA_END()
#endif

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
                             const boost::shared_ptr<Constant>& C)
    : RayleighImpBase(Pres, Alt, C)
{
}

//-----------------------------------------------------------------------
/// Calculate the rayleigh cross section for the given
/// wavenumber/wavelength.
//-----------------------------------------------------------------------

DoubleWithUnit RayleighYoung::cross_section(const DoubleWithUnit& W) const
{
    double wn = W.convert_wave(units::inv_cm).value;

    const double a = constants->rayleigh_a().value;
    const double b = constants->rayleigh_b().value;
    const double depolar_fact = constants->rayleigh_depolarization_factor();

    double rayleigh_cross_section;
    rayleigh_wrap(&wn, &depolar_fact, &a, &b, &rayleigh_cross_section);

    return DoubleWithUnit(rayleigh_cross_section, Unit("m^2"));
}

//-----------------------------------------------------------------------
/// Make a copy of this class
//-----------------------------------------------------------------------

boost::shared_ptr<Rayleigh> RayleighYoung::clone() const
{
    return boost::shared_ptr<Rayleigh>(new RayleighYoung(pres, alt, constants));
}
