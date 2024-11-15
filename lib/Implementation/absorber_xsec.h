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

    // Density values
    ArrayAdWithUnit<double, 1> air_density;
    ArrayAdWithUnit<double, 2> gas_density;

private:
    void cache_height_delta(const boost::shared_ptr<Pressure>& press, const std::vector<boost::shared_ptr<Altitude> >& alt); 
    void cache_total_air_number_density_level();
    void cache_gas_number_density_level(const boost::shared_ptr<Pressure>& press, const std::vector<boost::shared_ptr<AbsorberVmr> >& vmr);

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
    /// Number density of air at each level boundary 
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 1> total_air_number_density_level() const;

    //-----------------------------------------------------------------------
    /// Number density of gas molecules at each level boundary
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 2> gas_number_density_level() const;

    //-----------------------------------------------------------------------
    // See base class description 
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 1> total_air_number_density_layer(int spec_index) const;

    //-----------------------------------------------------------------------
    // See base class description 
    //-----------------------------------------------------------------------

    virtual ArrayAdWithUnit<double, 2> gas_number_density_layer(int spec_index) const;

    //-----------------------------------------------------------------------
    // See base class description 
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 2> optical_depth_each_layer(double wn, int spec_index) const;

    virtual void print(std::ostream& Os) const;
    virtual boost::shared_ptr<Absorber> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& Gas_name) const;
    virtual boost::shared_ptr<XSecTable> xsec_table(const std::string& Gas_name) const;

    const boost::shared_ptr<Pressure> pressure() const
    {
        return press;
    }

    const boost::shared_ptr<Temperature> temperature() const
    {
        return temp;
    }

protected:
  // Directors with swig seems to need to have this available
  AbsorberXSec() = default;
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

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
};

} // namespace FullPhysics  

FP_EXPORT_KEY(AbsorberXSecCache);
FP_EXPORT_KEY(AbsorberXSec);
#endif
