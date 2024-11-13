#ifndef RT_DRIVER_SCENARIO_H
#define RT_DRIVER_SCENARIO_H

#include <blitz/array.h>

#include "printable.h"
#include "array_ad.h"

namespace FullPhysics {

/****************************************************************//**
 *******************************************************************/

class RtDriverScenarioGeometry : public Printable<RtDriverScenarioGeometry> {

public:

    RtDriverScenarioGeometry(double sza, double zen, double azm);

    double solar_zenith() { return sza_; }
    void solar_zenith(double sza) { sza_ = sza; }

    double observation_zenith() { return zen_; }
    void observation_zenith(double zen) { zen_ = zen; }

    double relative_azimuth() { return azm_; }
    void relative_azimuth(double azm) { azm_ = azm; }

    void check_configuration() const;

    virtual void print(std::ostream& Os, bool Short_form = false) const;
  
private:

    double sza_;
    double zen_;
    double azm_;

};

/****************************************************************//**
 *******************************************************************/

class RtDriverScenarioSurface : public Printable<RtDriverScenarioSurface> {

public:

    RtDriverScenarioSurface(int surface_type, bool do_thermal);

    int surface_type() const { return surface_type_; }
    void surface_type(const int surf_type) { surface_type_ = surf_type; }

    const bool has_thermal() const { return do_thermal_; }

    const int number_jacobian_parameter() const { return jac_perturbation_.rows(); }

    const ArrayAd<double, 1> surface_parameters() const { return surface_params_; }
    void surface_parameters(const ArrayAd<double, 1>& surface_params) { surface_params_ = surface_params; }

    double black_body() const { return black_body_; }
    void black_body(const double bb) { black_body_ = bb; }

    const blitz::Array<double, 1> jacobian_perturbation() { return jac_perturbation_; }
    void jacobian_perturbation(const blitz::Array<double, 1>& perturbation) { jac_perturbation_ = perturbation; }

    // Initalize simple example values
    void compute_simple_black_body(double wn = 568.69, double temperature = 290.0);

    void check_configuration() const;

    virtual void print(std::ostream& Os, bool Short_form = false) const;
  
private:

    int surface_type_;
    ArrayAd<double, 1> surface_params_;

    double black_body_;
    bool do_thermal_;

    blitz::Array<double, 1> jac_perturbation_;

};

/****************************************************************//**
 *******************************************************************/

class RtDriverScenarioAtmosphere : public Printable<RtDriverScenarioAtmosphere> {

public:

    RtDriverScenarioAtmosphere(int nlayer, int nmoms, int nstokes, bool do_aerosol, bool do_thermal);
  
    const bool has_aerosol() const { return pf_.rows() > 0; }
    const bool has_thermal() const { return black_body_.rows() > 0; }

    const int number_layer() const { return taug_.rows(); }
    const int number_moments() const { return pf_.rows() - 1; }
    const int number_stokes() const { return pf_.depth(); }
    const int number_jacobian_parameter() const { return jac_perturbation_.rows(); }
  
    // Get/set methods, done this way mainly for Python interface
    const blitz::Array<double, 1> height_grid() { return heights_; }
    void height_grid(const blitz::Array<double, 1>& heights) { heights_ = heights; }

    const ArrayAd<double, 1> od_gas() const { return taug_; }
    void od_gas(const ArrayAd<double, 1>& taug) { taug_ = taug; }

    const ArrayAd<double, 1> od_rayleigh() const { return taur_; }
    void od_rayleigh(const ArrayAd<double, 1>& taur) { taur_ = taur; }

    const ArrayAd<double, 1> od_aerosol() const { return taua_; }
    void od_aerosol(const ArrayAd<double, 1>& taua) { taua_ = taua; }

    const ArrayAd<double, 1> total_od() const;
    const ArrayAd<double, 1> singe_scattering_albedo(double aer_prop_ssa) const;

    const ArrayAd<double, 3> phase_function() const { return pf_; }
    void phase_function(const ArrayAd<double, 3>& pf) { pf_ = pf; }

    const blitz::Array<double, 1> black_body() { return black_body_; }
    void black_body(const blitz::Array<double, 1>& bb) { black_body_ = bb; }

    /// Perturbations for finite difference jacobians
    const blitz::Array<double, 1> jacobian_perturbation() { return jac_perturbation_; }
    void jacobian_perturbation(const blitz::Array<double, 1>& perturbation) { jac_perturbation_ = perturbation; }

    // Initialize simple example values
    void compute_simple_pf(double aer_prop_ssa = 0.9, double aer_prop_asym = 0.7, double depol = 0.0);
    void compute_simple_black_body(double wn = 568.69, double temperature = 290.0);
    void compute_simple_height_grid();

    // Set jacobian values for taug, taur, taua
    void initialize_jacobians();

    void check_configuration() const;
  
    virtual void print(std::ostream& Os, bool Short_form = false) const;
    
private:
  

    blitz::Array<double, 1> heights_;

    ArrayAd<double, 1> taug_;
    ArrayAd<double, 1> taur_;
    ArrayAd<double, 1> taua_;
    ArrayAd<double, 3> pf_;

    blitz::Array<double, 1> black_body_;
  
    blitz::Array<double, 1> jac_perturbation_;

};

}

#endif
