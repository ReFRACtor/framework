#ifndef PCA_OPTICAL_PROP_H
#define PCA_OPTICAL_PROP_H

#include "generic_object.h"
#include "atmosphere_standard.h"
#include "ground_lambertian.h"
#include "spectral_domain.h"
#include "hdf_file.h"

namespace FullPhysics {

/****************************************************************//**
  Represents the interface of the values needed by PCA to perform
  its calculations. The values here are for multiple spectral points
  across a whole channel/band.
 *******************************************************************/
class PCAOpticalProperties : public virtual GenericObject {
public:
    virtual ~PCAOpticalProperties() = default;

    virtual blitz::Array<double, 2> gas_optical_depth() const = 0;
    virtual blitz::Array<double, 2> total_optical_depth() const = 0;
    virtual blitz::Array<double, 2> single_scattering_albedo() const = 0;
    virtual blitz::Array<int, 1> primary_gas_dominates() const = 0;
    virtual blitz::Array<double, 3> intermediate_variable() const = 0;
    virtual blitz::Array<double, 1> surface_albedo() const = 0;

};

/****************************************************************//**
  Computes atmospheric optical properties need repeatedly in PCA
  calculations. Uses the framework's atmosphere class to gather
  the values.
 *******************************************************************/

class PCAOpticalPropertiesAtmosphere : public PCAOpticalProperties {
public:

    PCAOpticalPropertiesAtmosphere(const boost::shared_ptr<AtmosphereStandard>& atm, const SpectralDomain& spec_domain, int channel_index, std::string primary_absorber, bool show_progress=true);

    virtual blitz::Array<double, 1> wavenumber() const { return wavenumber_; }
    virtual blitz::Array<double, 2> gas_optical_depth() const { return gas_optical_depth_; }
    virtual blitz::Array<double, 2> total_optical_depth() const { return total_optical_depth_; }
    virtual blitz::Array<double, 2> single_scattering_albedo() const { return single_scattering_albedo_; }
    virtual blitz::Array<int, 1> primary_gas_dominates() const { return primary_gas_dominates_; }
    virtual blitz::Array<double, 3> intermediate_variable() const { return intermediate_; }
    virtual blitz::Array<double, 1> surface_albedo() const { return surface_albedo_; }

private:

    void compute_properties();

    boost::shared_ptr<AtmosphereStandard> atmosphere;
    boost::shared_ptr<GroundLambertian> lambertian;

    bool show_progress_;

    int channel_index_;
    int primary_abs_index_;

    blitz::Array<double, 1> wavenumber_;
    blitz::Array<double, 2> gas_optical_depth_;
    blitz::Array<double, 2> total_optical_depth_;
    blitz::Array<double, 2> single_scattering_albedo_;
    blitz::Array<int, 1> primary_gas_dominates_;
    blitz::Array<double, 3> intermediate_;
    blitz::Array<double, 1> surface_albedo_;
};

/****************************************************************//**
  Loads precomputed atmospheric optical properties from a file.
  Mainly intended for use in testing.
 *******************************************************************/


class PCAOpticalPropertiesFile : public PCAOpticalProperties {
public:
    PCAOpticalPropertiesFile(const std::string& filename);
    PCAOpticalPropertiesFile(const boost::shared_ptr<HdfFile>& file);

    virtual blitz::Array<double, 2> gas_optical_depth() const;
    virtual blitz::Array<double, 2> total_optical_depth() const;
    virtual blitz::Array<double, 2> single_scattering_albedo() const;
    virtual blitz::Array<int, 1> primary_gas_dominates() const;
    virtual blitz::Array<double, 3> intermediate_variable() const;
    virtual blitz::Array<double, 1> surface_albedo() const;

private:
    boost::shared_ptr<HdfFile> hdf_file;

};

}

#endif
