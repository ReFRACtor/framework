#include "rayleigh_bodhaine.h"
#include "fp_serialize_support.h"

#include <cmath>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void RayleighBodhaine::serialize(Archive & ar,
				 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RayleighImpBase);
}

FP_IMPLEMENT(RayleighBodhaine);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(RayleighBodhaine, Rayleigh)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
                          const std::vector<boost::shared_ptr<Altitude> >&,
                          const boost::shared_ptr<Constant>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
RayleighBodhaine::RayleighBodhaine(const boost::shared_ptr<Pressure>& Pres,
                                   const std::vector<boost::shared_ptr<Altitude> >& Alt,
                                   const boost::shared_ptr<Constant>& C)
    : RayleighImpBase(Pres, Alt, C)
{
}

//-----------------------------------------------------------------------
/// Calculate the rayleigh cross section for the given
/// wavenumber/wavelength.
//-----------------------------------------------------------------------

DoubleWithUnit RayleighBodhaine::cross_section(const DoubleWithUnit& W) const
{
    double wl_um = W.convert_wave(units::micron).value;
    double wl_nm = W.convert_wave(units::nm).value;

    double wl_cm = wl_nm * 1.0e-7; // wl [cm]
    double wl_um_m2 = 1.0 / std::pow(wl_um, 2.0);

    // Avogadro's number [molecules/mol]
    const double A = 6.0221367e23;

    // molar volume at 273.15 K and 1013.25 mbar [L/mol]
    const double mvol = 22.4141e0;

    // number density (at 288.15 K and 1013.25 mbar)
    // number density [molecules/cm^3]
    const double numdens = A * 273.15 / (mvol * 288.15) / 1.0e3;

    // vmr
    const double pvmr[4] = {78.084, 20.946, 0.934, 0.036};
    const double pvmr_sum = pvmr[0] + pvmr[1] + pvmr[2] + pvmr[3];

    // King factor for depolarization
    // Earth (Bodhaine et al. [1999])
    // King factor = (6+3*rho)/(6-7*rho)
    Array<double, 1> F(4);
    F(0) = 1.034 + 3.17e-4 * wl_um_m2;
    F(1) = 1.096 + 1.385e-3 * wl_um_m2 + 1.448e-4 * std::pow(wl_um_m2, 2.0);
    F(2) = 1.0;
    F(3) = 1.15;

    // sum(pvmr(:)*F(:))/sum(pvmr)
    double king_factor = 
        (pvmr[0] * F(0) + pvmr[1] * F(1) + pvmr[2] * F(2) + pvmr[3] * F(3)) / pvmr_sum;

    // refractive index (at 288.15 K and 1013.25 mbar)
    // Earth (Bodhaine et al. [1999])
    double refr = 1.0 + (8060.51 + 2480990.0 / (132.274 - wl_um_m2) +
            17455.7 / (39.32957 - wl_um_m2)) *1.0e-8 * 
            (1.0 + 0.54 * (pvmr[3] / 100.0 -0.0003));

    double index_factor = std::pow((std::pow(refr, 2.0) - 1.0) / (std::pow(refr, 2.0) + 2.0), 2.0);

    // cross section in cm^2 per molecule of gas
    double ray_ext = 24.0 * std::pow(M_PI, 3.0) / std::pow(numdens, 2.0) / std::pow(wl_cm, 4.0) * index_factor * king_factor;

    return DoubleWithUnit(ray_ext, Unit("cm^2"));
}

//-----------------------------------------------------------------------
/// Make a copy of this class
//-----------------------------------------------------------------------

boost::shared_ptr<Rayleigh> RayleighBodhaine::clone() const
{
    return boost::shared_ptr<Rayleigh>(new RayleighBodhaine(pres, alt, constants));
}
