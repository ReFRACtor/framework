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

class AbsorberXSec;
class AbsorberXSecCache : public CalculationCache<AbsorberXSec> {
public:
    AbsorberXSecCache() = default;
    ~AbsorberXSecCache() = default;

    virtual void fill_cache(const AbsorberXSec& absorber);

    // pressure / temperature grid computed from Pressure / Temperature objects
    ArrayAdWithUnit<double, 1> pgrid;
    ArrayAdWithUnit<double, 1> tgrid;

    // Height of each layer, the difference in the level altitude values
    std::vector<ArrayAdWithUnit<double, 1> > height_delta_layer;

private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  This class maintains the absorber portion of the state. This
  particular implementation uses UV Cross Section (XSec) tables.
*******************************************************************/

class AbsorberXSec: virtual public Absorber {
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
    /// Wet air density at each level boundary determined using Loschmidt's constant
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 1> air_density_level() const;

    //-----------------------------------------------------------------------
    /// Number density of the the gas molecule at each level boundary
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

    // Cached values that only get recomputed when the cache is invalidated
    friend AbsorberXSecCache;
    mutable AbsorberXSecCache cache;

    AbsorberXSec() = default;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
};

} // namespace FullPhysics  

FP_EXPORT_KEY(AbsorberXSecCache);
FP_EXPORT_KEY(AbsorberXSec);
#endif
