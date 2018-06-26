#ifndef GAS_VMR_APRIORI_H
#define GAS_VMR_APRIORI_H

#include <blitz/array.h>

#include "meteorology.h"
#include "level_1b_sample_coefficient.h"
#include "altitude.h"
#include "hdf_file.h"
#include "pressure.h"

#include "reference_vmr_apriori.h"

namespace FullPhysics {

/****************************************************************//**
 Adapts the ReferenceVmrApriori class into a form that is 
 easier to work with in the context of how this framework works.

 This class deals keys off of pressure levels like the rest of
 the framework. It also uses the increasing pressure levels 
 convention.

 The VMR returned will be in the order of levels expected 
 elsewhere in the framework.
*******************************************************************/

class GasVmrApriori : public Printable<GasVmrApriori> {
public:
    GasVmrApriori(const boost::shared_ptr<Meteorology>& met_file,
                  const boost::shared_ptr<Level1bSampleCoefficient>& l1b_file,
                  const boost::shared_ptr<Altitude>& alt,
                  const HdfFile& hdf_static_input,
                  const std::string& hdf_group,
                  const std::string& gas_name,
                  const int temp_avg_window = 11);

    GasVmrApriori(const blitz::Array<double, 1>& pressure_levels,
                  const blitz::Array<double, 1>& temperature_levels,
                  const double& obs_latitude,
                  const Time& obs_time,
                  const boost::shared_ptr<Altitude>& altitude,
                  const HdfFile& reference_file,
                  const std::string& hdf_group,
                  const std::string& gas_name,
                  const int temp_avg_window = 11);

    const blitz::Array<double, 1> apriori_vmr() const;
    const blitz::Array<double, 1> apriori_vmr(const Pressure& pressure) const;

    const boost::shared_ptr<ReferenceVmrApriori> reference() const { return ref_apriori; }

    const DoubleWithUnit tropopause_altitude() const { return ref_apriori->model_tropopause_altitude(); }
    const double tropopause_pressure() const;

    void print(std::ostream& Os) const { Os << "GasVmrApriori"; }

private:

    void initialize(const blitz::Array<double, 1>& pressure_levels,
                    const blitz::Array<double, 1>& temperature_levels,
                    const double& obs_latitude,
                    const Time& obs_time,
                    const boost::shared_ptr<Altitude>& altitude,
                    const HdfFile& reference_file,
                    const std::string& hdf_group,
                    const std::string& gas_name,
                    const int temp_avg_window);

    boost::shared_ptr<ReferenceVmrApriori> ref_apriori;

    blitz::Array<double, 1> model_alt;
    blitz::Array<double, 1> model_press;
    blitz::Array<double, 1> ref_vmr;

    std::string gas_name;
};

}

#endif
