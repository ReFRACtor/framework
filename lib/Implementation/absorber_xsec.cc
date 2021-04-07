#include "fp_serialize_support.h"
#include "ostream_pad.h"

#include "absorber_xsec.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberXSecCache::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CacheInvalidatedObserver);
}

template<class Archive>
void AbsorberXSec::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Absorber)
    & FP_NVP(press) & FP_NVP(temp) & FP_NVP(alt) & FP_NVP(vmr);
}

FP_IMPLEMENT(AbsorberXSecCache);
FP_IMPLEMENT(AbsorberXSec);
#endif

AbsorberXSec::AbsorberXSec(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
                           const boost::shared_ptr<Pressure>& Press,
                           const boost::shared_ptr<Temperature>& Temp,
                           const std::vector<boost::shared_ptr<Altitude> >& Alt,
                           const std::vector<boost::shared_ptr<XSecTable> >& XSec_tables)
: press(Press), temp(Temp), alt(Alt), vmr(Vmr), xsec_tables(XSec_tables)
{
    press->add_cache_invalidated_observer(cache);
    temp->add_cache_invalidated_observer(cache);
}

void AbsorberXSecCache::fill_cache(const AbsorberXSec& absorber)
{
    pgrid.reference(absorber.press->pressure_grid().convert(units::Pa));
    tgrid.reference(absorber.temp->temperature_grid(*(absorber.press)));

    // Number of derivative variables for auto derivatives
    int num_jac_variable = pgrid.number_variable();

    height_delta_layer.clear();

    for(int sensor_idx = 0; sensor_idx < absorber.number_spectrometer(); sensor_idx++) {
        ArrayAd<double, 1> height_delta_sensor(absorber.press->number_layer(), num_jac_variable);

        Unit height_units;
        for(int lay_idx = 0; lay_idx < absorber.press->number_layer(); lay_idx++) {
            AutoDerivativeWithUnit<double> h_lay_top = absorber.alt[sensor_idx]->altitude(pgrid(lay_idx));
            AutoDerivativeWithUnit<double> h_lay_bottom = absorber.alt[sensor_idx]->altitude(pgrid(lay_idx+1));
            height_delta_sensor(lay_idx) = h_lay_top.value - h_lay_bottom.value;
            height_units = h_lay_top.units;
        }

        height_delta_layer.push_back(ArrayAdWithUnit<double, 1>(height_delta_sensor, height_units));
    }

}

ArrayAdWithUnit<double, 1> AbsorberXSec::air_density_level() const
{
    cache.fill_cache_if_needed(*this);

    // Loschmidt's number (particles/cm3), STP parameters
    const DoubleWithUnit rho_stand(2.68675e+19, "cm^-3");
    const DoubleWithUnit pzero(1013.25e2, units::Pa); // Standard presure
    const DoubleWithUnit tzero(273.15e0, units::K);   // Standard temperature 0 degC
    const DoubleWithUnit rho_zero = rho_stand * tzero / pzero;
    const DoubleWithUnit dens_const = 1.0e+05 * rho_zero;

    ArrayAd<double, 1> air_density(press->number_level(), cache.pgrid.number_variable());;
    for(int lev_idx = 0; lev_idx < press->number_level(); lev_idx++) {
        air_density(lev_idx) = dens_const.value * cache.pgrid.value(lev_idx) / cache.tgrid.value(lev_idx);
    }

    // Units here will be cm^-3 since the press/temp ratio cancels out with that in the constant
    // Set the unit explicitly so there is not "cruft" factors left over from using the Unit class
    return ArrayAdWithUnit<double, 1>(air_density, "cm^-3");
}

ArrayAdWithUnit<double, 2> AbsorberXSec::gas_density_level() const
{
    cache.fill_cache_if_needed(*this);

    ArrayAdWithUnit<double, 1> air_density(air_density_level());

    ArrayAd<double, 2> gas_density(air_density.rows(), vmr.size(), air_density.number_variable());

    for(int gas_idx = 0; gas_idx < vmr.size(); gas_idx++) {
        ArrayAd<double, 1> vmr_grid(vmr[gas_idx]->vmr_grid(*press));

        if(gas_density.number_variable() == 0 and vmr_grid.number_variable() > 0) {
            gas_density.resize_number_variable(vmr_grid.number_variable());
        }

        for(int lev_idx = 0; lev_idx < press->number_level(); lev_idx++) {
            gas_density(lev_idx, gas_idx) = air_density.value(lev_idx) * vmr_grid(lev_idx);
        }
    }

    return ArrayAdWithUnit<double, 2>(gas_density, air_density.units);
}

ArrayAd<double, 2> AbsorberXSec::optical_depth_each_layer(double wn, int spec_index) const
{
    cache.fill_cache_if_needed(*this);

    range_check(spec_index, 0, number_spectrometer());

    ArrayAdWithUnit<double, 2> gas_density(gas_density_level());
    DoubleWithUnit spectral_point(wn, units::inv_cm);

    ArrayAd<double, 2> gas_od(press->number_layer(), vmr.size(), gas_density.number_variable());

    for(int gas_idx = 0; gas_idx < vmr.size(); gas_idx++) {
        ArrayAd<double, 1> od_unweighted(
            xsec_tables[gas_idx]->optical_depth_each_layer_unweighted(spectral_point, gas_density.value(Range::all(), gas_idx), cache.tgrid.value)
        );

        for(int lay_idx = 0; lay_idx < press->number_layer(); lay_idx++) {
            AutoDerivativeWithUnit<double> height_diff(cache.height_delta_layer[spec_index](lay_idx));
            gas_od(lay_idx, gas_idx) = height_diff.value * od_unweighted(lay_idx); 
        }
    }

    return gas_od;
}

void AbsorberXSec::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "AbsorberXSec:\n";
    for(int i = 0; i < (int) vmr.size(); ++i) {
        Os << "  Absorber VMR[" << i << "]:\n";
        opad << *vmr[i] << "\n";
        opad.strict_sync();
    }
}

boost::shared_ptr<Absorber> AbsorberXSec::clone() const
{
    return boost::shared_ptr<Absorber>(new AbsorberXSec(vmr, press, temp, alt, xsec_tables));
}

boost::shared_ptr<AbsorberVmr> AbsorberXSec::absorber_vmr(const std::string& Gas_name) const
{
    int i = gas_index(Gas_name);
    if(i < 0) {
        Exception err;
        err << "Gas named " << Gas_name << " not present in vmr list";
        throw err;
    }
    return vmr[i];
}
