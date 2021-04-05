#ifndef ABSORBER_XSEC_H
#define ABSORBER_XSEC_H
#include <vector>

#include "fp_exception.h"

#include "absorber.h"
#include "absorber_vmr.h"
#include "pressure.h"
#include "temperature.h"
#include "altitude.h"
#include "constant.h"
#include "calculation_cache.h"

#include "xsec_table.h"

namespace FullPhysics {

/****************************************************************//**
  This class maintains the absorber portion of the state. This
  particular implementation uses UV Cross Section (XSec) tables.
*******************************************************************/

class AbsorberXSec: virtual public Absorber,
                            public Observer<AbsorberVmr>,
                            public Observer<Pressure>,
                            public Observer<Temperature>,
                            public Observer<Altitude>,
                            public CalculationCache<AbsorberXSec> {
public:
    AbsorberXSec(const std::vector<boost::shared_ptr<AbsorberVmr> > Vmr,
                 const boost::shared_ptr<Pressure>& Press,
                 const boost::shared_ptr<Temperature>& Temp,
                 const std::vector<boost::shared_ptr<Altitude> >& Alt,
                 const std::vector<boost::shared_ptr<XSecTable> >& XSec_tables);

    virtual ~AbsorberXSec() = default;

    virtual int number_species() const
    {
        return (int) vmr.size();
    }
    virtual int number_spectrometer() const
    {
        return (int) alt.size();
    }
    virtual int number_layer() const
    {
        return press->number_layer();
    }
    virtual std::string gas_name(int Species_index) const
    {
        range_check(Species_index, 0, number_species());
        return vmr[Species_index]->gas_name();
    }

    //-----------------------------------------------------------------------
    /// For performance, we cache some data as we calculate it. This
    /// becomes stale when the pressure is changed, so we observe press
    /// and mark the cache when it changes.
    //-----------------------------------------------------------------------

    virtual void notify_update(const Pressure& UNUSED(P) )
    {
        invalidate_cache();
    }
    virtual void notify_update(const Temperature& UNUSED(T))
    {
        invalidate_cache();
    }
    virtual void notify_update(const Altitude& UNUSED(A))
    {
        invalidate_cache();
    }
    virtual void notify_update(const AbsorberVmr& UNUSED(A) )
    {
        invalidate_cache();
    }

    virtual void fill_cache(const AbsorberXSec& T);

    //-----------------------------------------------------------------------
    /// Dry air density on each level determined using Loschmidt's constant
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 1> air_density_level() const;

    //-----------------------------------------------------------------------
    /// Number density of the the gas molecule per level
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 2> gas_density_level() const;

    //-----------------------------------------------------------------------
    // See base class description 
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 2> optical_depth_each_layer(double wn, int spec_index) const;

    virtual void print(std::ostream& Os) const;
    virtual boost::shared_ptr<Absorber> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& Gas_name) const;

    const Pressure& pressure() const
    {
        return *press;
    }

private:
    // Objects used to calculate the optical depth
    boost::shared_ptr<Pressure> press;
    boost::shared_ptr<Temperature> temp;
    std::vector<boost::shared_ptr<Altitude> > alt;
    std::vector<boost::shared_ptr<AbsorberVmr> > vmr;

    // Filenames contain cross section data
    std::vector<boost::shared_ptr<XSecTable> > xsec_tables;

    // Scratch variable used to calculate taug. We keep this around so
    // we don't keep recreating the Array (so this is to improve performance)
    mutable ArrayAd<double, 2> taug;

    AbsorberXSec() = default;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
};

} // namespace FullPhysics  

FP_EXPORT_KEY(AbsorberXSec);
#endif
