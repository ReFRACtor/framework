#include "rayleigh_imp_base.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

RayleighImpBase::RayleighImpBase(const boost::shared_ptr<Pressure>& Pres,
                                 const std::vector<boost::shared_ptr<Altitude> >& Alt,
                                 const boost::shared_ptr<Constant>& C)
    : constants(C), pres(Pres), alt(Alt), cache_is_stale(true)
{
    pres->add_observer(*this);
    BOOST_FOREACH(boost::shared_ptr<Altitude>& a, alt) {
        a->add_observer(*this);
    }
}

//-----------------------------------------------------------------------
/// This gives the optical depth for each layer, for the given wave
/// number. Note this only includes the Rayleigh portion of this,
/// Atmosphere class combines this with Absorbers and Aerosol
/// scattering.
///
/// This has size of pres->number_active_layer().
//-----------------------------------------------------------------------

ArrayAd<double, 1> RayleighImpBase::optical_depth_each_layer(double wn, int spec_index) const
{
    range_check(spec_index, 0, (int) alt.size());
    fill_cache();

    ArrayAd<double, 1> res(part_independent_wn(spec_index, Range::all()).copy());
    double rayleigh_cross_section = cross_section(DoubleWithUnit(wn, units::inv_cm)).convert(Unit("m^2")).value;

    res.value() *= rayleigh_cross_section;
    res.jacobian() *= rayleigh_cross_section;
    return res;
}

//-----------------------------------------------------------------------
/// Fill in cache, if needed.
//-----------------------------------------------------------------------
void RayleighImpBase::fill_cache() const
{
    int local_nvar = part_independent_wn.number_variable();

    int alt_nvar = 0;

    for(int i = 0; i < (int) alt.size(); ++i) {
        alt_nvar = std::max(alt_nvar, alt[i]->gravity(pres->pressure_grid()(0)).value.number_variable());
    }

    if(!cache_is_stale && local_nvar == alt_nvar) {
        return;
    }

    // A comment in the old fortran code indicates this is Avogadro's
    // number. However, this actually is not (Avogadro's number is
    // 6.02214179e23 /mol.  I'm not actually sure what this number is.
    // You can look at the file exe/full_physics/src/rtmod/ray_tau.F90
    // in tagged B2.06.02_plus_2.07_backport version to find the
    // original code and comment
    const double a0 = 6.02297e26;

    const double molar_weight_dry_air = constants->molar_weight_dry_air().convert(Unit("g / mol")).value;

    int nvar = std::max(pres->pressure_grid().value.number_variable(), alt_nvar);
    part_independent_wn.resize((int) alt.size(), pres->number_layer(), nvar);

    for(int i = 0; i < part_independent_wn.cols(); ++i)
        for(int j = 0; j < part_independent_wn.rows(); ++j) {
            AutoDerivativeWithUnit<double> deltap = pres->pressure_grid()(i + 1) -
                                                    pres->pressure_grid()(i);
            AutoDerivativeWithUnit<double> play =
                pres->pressure_grid()(i) + deltap / 2;
            part_independent_wn(j, i) = a0 * deltap.convert(units::Pa).value /
                                        (molar_weight_dry_air * alt[j]->gravity(play).convert(Unit("m/s^2")).value);
        }

    cache_is_stale = false;
}
