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
    & FP_NVP(press) & FP_NVP(temp) & FP_NVP(alt) & FP_NVP(vmr) & FP_NVP(xsec_tables);
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
    if(vmr.size() == 0) {
        throw Exception("No vmr objects passed");
    }

    if(xsec_tables.size() == 0) {
        throw Exception("No cross section tables passed");
    }

    if(vmr.size() != xsec_tables.size()) {
        Exception err;
        err << "Number of AbsorberVmr objects: " << vmr.size()
            << " does not match the number of XSecTable objects: " << xsec_tables.size();
        throw err;
    }

    // Ensure that items involved in cache computation trigger a recomputation if updated
    press->add_cache_invalidated_observer(cache);
    temp->add_cache_invalidated_observer(cache);

    for(auto& alti : alt)
      alti->add_cache_invalidated_observer(cache);

}

void AbsorberXSecCache::fill_cache(const AbsorberXSec& absorber)
{
    pgrid.reference(absorber.press->pressure_grid().convert(units::Pa));
    tgrid.reference(absorber.temp->temperature_grid(*(absorber.press)).convert(units::K));

    cache_height_delta(absorber.press, absorber.alt);
    cache_total_air_number_density_level();
    cache_gas_number_density_level(absorber.press, absorber.vmr);
}

void AbsorberXSecCache::cache_height_delta(const boost::shared_ptr<Pressure>& press, const std::vector<boost::shared_ptr<Altitude> >& alt)
{
    // Number of derivative variables for auto derivatives
    int num_jac_variable = pgrid.number_variable();

    height_delta_layer.clear();

    for(auto& alti : alt) {
        ArrayAd<double, 1> height_delta_sensor(press->number_layer(), num_jac_variable);

        Unit height_units;
        for(int lay_idx = 0; lay_idx < press->number_layer(); lay_idx++) {
            AutoDerivativeWithUnit<double> h_lay_top = alti->altitude(pgrid(lay_idx));
            AutoDerivativeWithUnit<double> h_lay_bottom = alti->altitude(pgrid(lay_idx+1));
            height_delta_sensor(lay_idx) = h_lay_top.value - h_lay_bottom.value;
            height_units = h_lay_top.units;
        }

        height_delta_layer.push_back(ArrayAdWithUnit<double, 1>(height_delta_sensor, height_units));
    }
}

void AbsorberXSecCache::cache_total_air_number_density_level()
{
    // Loschmidt's number (particles/cm3), STP parameters
    const DoubleWithUnit rho_stand(2.68675e+19, "cm^-3");
    const DoubleWithUnit pzero(1013.25e2, units::Pa); // Standard presure
    const DoubleWithUnit tzero(273.15e0, units::K);   // Standard temperature 0 degC
    const DoubleWithUnit rho_zero = rho_stand * tzero / pzero;

    if (pgrid.rows() == 0) {
        throw Exception("pgrid not yet computed by cache");
    }

    air_density.value.resize(pgrid.rows(), pgrid.number_variable());;

    // Units here will be cm^-3 since the press/temp ratio cancels out with that in the constant
    air_density.units = Unit("cm^-3");

    for(int lev_idx = 0; lev_idx < pgrid.rows(); lev_idx++) {
        air_density.value(lev_idx) = rho_zero.value * pgrid.value(lev_idx) / tgrid.value(lev_idx);
    }
}

void AbsorberXSecCache::cache_gas_number_density_level(const boost::shared_ptr<Pressure>& press, const std::vector<boost::shared_ptr<AbsorberVmr> >& vmr)
{
    if(air_density.rows() == 0) {
        throw Exception("air_density not yet computed by cache");
    }

    gas_density.value.resize(air_density.rows(), vmr.size(), air_density.number_variable());
    gas_density.units = air_density.units;

    for(int gas_idx = 0; gas_idx < gas_density.value.cols(); ++gas_idx) {
      ArrayAd<double, 1> vmr_grid(vmr[gas_idx]->vmr_grid(*press));
      if(gas_density.number_variable() == 0 and vmr_grid.number_variable() > 0) {
	gas_density.value.resize_number_variable(vmr_grid.number_variable());
      }

      for(int lev_idx = 0; lev_idx < press->number_level(); lev_idx++) {
	gas_density.value(lev_idx, gas_idx) = air_density.value(lev_idx) * vmr_grid(lev_idx);
      }
    }
}

ArrayAdWithUnit<double, 1> AbsorberXSec::total_air_number_density_level() const
{
    cache.fill_cache_if_needed(*this);
    return cache.air_density;
}

ArrayAdWithUnit<double, 2> AbsorberXSec::gas_number_density_level() const
{
    cache.fill_cache_if_needed(*this);
    return cache.gas_density;
}

ArrayAdWithUnit<double, 1> AbsorberXSec::total_air_number_density_layer(int spec_index) const
{
    cache.fill_cache_if_needed(*this);

    // Ensure consistency of units of items to be integraded
    ArrayAdWithUnit<double, 1> height_delta(cache.height_delta_layer[spec_index].convert(Unit("cm")));
    ArrayAdWithUnit<double, 1> air_dens_lev(total_air_number_density_level().convert(Unit("cm^-3")));

    ArrayAd<double, 1> air_dens_lay(press->number_layer(), air_dens_lev.number_variable());
   
    for(int lay_idx = 0; lay_idx < press->number_layer(); lay_idx++) {
        air_dens_lay(lay_idx) = height_delta(lay_idx).value * (air_dens_lev(lay_idx).value + air_dens_lev(lay_idx+1).value) * 0.5;
    }

    // One cm^-1 cancels out with altitude units
    return ArrayAdWithUnit<double, 1>(air_dens_lay, Unit("cm^-2"));
}

ArrayAdWithUnit<double, 2> AbsorberXSec::gas_number_density_layer(int spec_index) const
{
    cache.fill_cache_if_needed(*this);

    // Ensure consistency of units of items to be integraded
    ArrayAdWithUnit<double, 1> height_delta(cache.height_delta_layer[spec_index].convert(Unit("cm")));
    ArrayAdWithUnit<double, 2> gas_dens_lev(gas_number_density_level().convert(Unit("cm^-3")));

    ArrayAd<double, 2> gas_dens_lay(press->number_layer(), number_species(), gas_dens_lev.number_variable());
   
    for(int lay_idx = 0; lay_idx < press->number_layer(); lay_idx++)
      for(int gas_idx = 0; gas_idx < gas_dens_lay.cols(); gas_idx++)
	gas_dens_lay(lay_idx, gas_idx) = height_delta(lay_idx).value * (gas_dens_lev(lay_idx, gas_idx).value + gas_dens_lev(lay_idx+1, gas_idx).value) * 0.5;


    // One cm^-1 cancels out with altitude units
    return ArrayAdWithUnit<double, 2>(gas_dens_lay, Unit("cm^-2"));
}

ArrayAd<double, 2> AbsorberXSec::optical_depth_each_layer(double wn, int spec_index) const
{
    cache.fill_cache_if_needed(*this);

    range_check(spec_index, 0, number_spectrometer());

    DoubleWithUnit spectral_point(wn, units::inv_cm);

    ArrayAd<double, 2> gas_od(press->number_layer(), vmr.size(), cache.gas_density.number_variable());

    for(int gas_idx = 0; gas_idx < gas_od.cols(); gas_idx++) {
        ArrayAdWithUnit<double, 1> od_unweighted(
            xsec_tables[gas_idx]->optical_depth_each_layer_unweighted(spectral_point, cache.gas_density(Range::all(), gas_idx), cache.tgrid)
        );

        for(int lay_idx = 0; lay_idx < press->number_layer(); lay_idx++) {
            AutoDerivativeWithUnit<double> height_diff(cache.height_delta_layer[spec_index](lay_idx));
            gas_od(lay_idx, gas_idx) = height_diff.value * od_unweighted(lay_idx).convert(1/height_diff.units).value; 
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

boost::shared_ptr<XSecTable> AbsorberXSec::xsec_table(const std::string& Gas_name) const
{
    int i = gas_index(Gas_name);
    if(i < 0) {
        Exception err;
        err << "Gas named " << Gas_name << " not present in vmr list";
        throw err;
    }
    return xsec_tables[i];
}
