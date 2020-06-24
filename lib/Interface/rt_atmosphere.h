#ifndef RT_ATMOSPHERE_H
#define RT_ATMOSPHERE_H
#include "state_vector_observer.h"
#include "observer.h"
#include "array_ad_with_unit.h"
#include "accumulated_timer.h"
#include "ground.h"
#include "optical_properties.h"

namespace FullPhysics {
/****************************************************************//**
  This class is responsible for setting up the atmosphere and ground
  information needed to run the Radiative transfer code.

  There are many, many properties associated with the atmosphere. This
  class is not meant to model these properties, it is really the very
  limited information needed to run the Radiative transfer code.

  Note that this includes both the atmosphere and surface parameters
  needed by the RT code.

  The calculation of the Jacobians in LIDORT takes a time directly
  proportional to the number of variables we are taking the Jacobian
  with respect to, we use an "intermediate" set of variables for some
  of the reported gradients (e.g., AtmosphereStandard uses taur, taug, and
  tau for each of the aerosol). To support future Atmosphere classes, we
  are purposely vague on exactly what these intermediate variables are, at
  least through the RtAtmosphere interface. The OpticalProperties "intermediate_jacobian"
  function can be used to get the value of these intermediate
  variables and Jacobian with the state vector variables.

  A description of the intermediate variables can be found in
  doc/LIDORT_Jacobian.pdf.

  Other objects may depend on the RtAtmosphere, and should be updated
  when the RtAtmosphere is updated. To facilitate that, this class is
  an Oberverable, and objects can add themselves as Observers to be
  notified when the RtAtmosphere is updated.

  Because the absorber calculation tends to be a bottle neck, we keep
  a timer in this class. This class keeps track of the time used in
  the atmosphere calculations. Other classes can make use of
  this information for logging if desired.
*******************************************************************/

class RtAtmosphere : virtual public StateVectorObserver,
		     public Observable<RtAtmosphere> {
public:
  virtual ~RtAtmosphere() {}
  static AccumulatedTimer timer;
  virtual void add_observer(Observer<RtAtmosphere>& Obs)
  {
    add_observer_do(Obs, *this);
  }
  virtual void remove_observer(Observer<RtAtmosphere>& Obs)
  {
    remove_observer_do(Obs, *this);
  }
  virtual std::string timer_info() const;
  
  //-----------------------------------------------------------------------
  /// Number of layers we currently have.
  //-----------------------------------------------------------------------
  
  virtual int number_layer() const = 0;
  
  //-----------------------------------------------------------------------
  /// Number of spectrometers we have.
  //-----------------------------------------------------------------------
  
  virtual int number_spectrometer() const = 0;
  
  //-----------------------------------------------------------------------
  /// Altitude grid for current pressure grid.
  //-----------------------------------------------------------------------
  
  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index)
    const = 0;
  
  //-----------------------------------------------------------------------
  /// Total column optical depth for the given gas. This is 0 if the band
  /// isn't one that sees that gas.
  ///
  /// The jacobian output is with respect to the state vector.
  //-----------------------------------------------------------------------
  
  virtual AutoDerivative<double> column_optical_depth
  (double wn, int spec_index, const std::string& Gas_name) const = 0;

  //-----------------------------------------------------------------------
  /// The optical depth for each layer, for the given wave number.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// radiative transfer variables, rather than the state vector variables (see
  /// description of RtAtmosphere class for discussion of this).
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \return Optical depth for each layer. This is number_layer() in size
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> optical_depth_wrt_rt(double wn, int spec_index)
    const = 0;

  //-----------------------------------------------------------------------
  /// The single scattering albedo for each layer, for the given wave number.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// radiative transfer variables, rather than the state vector variables (see
  /// description of RtAtmosphere class for discussion of this).
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \return Single scattering albedo for each layer. This is
  /// number_layer() in size
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> single_scattering_albedo_wrt_rt
  (double wn, int spec_index) const = 0;

  //-----------------------------------------------------------------------
  /// The scattering moments for for each layer, for the given wave
  /// number.
  ///
  /// The scattering moments use the de Rooij convention for the 6
  /// scattering matrix element.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// radiative transfer variables, rather than the state vector variables (see
  /// description of RtAtmosphere class for discussion of this).
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \param nummom Number of moments to include in
  ///           scatt_mom_each_layer, the default it to include all of
  ///           them.
  /// \param numscat Number of scattering matrix elements to include in
  ///           scatt_mom_each_layer, the default it to include all of
  ///           them.
  /// \return Scattering moments for each layer. This is
  ///         number_moment + 1 x number_layer() x number scattering
  ///         matrix elements
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 3> phase_function_moments_wrt_rt
  (double wn, int spec_index, int nummom = -1, int numscat = -1) const = 0;

  //-----------------------------------------------------------------------
  /// The optical depth for each layer, for the given wave number.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// state vector variables.
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \return Optical depth for each layer. This is number_layer() in size
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> optical_depth_wrt_state_vector
  (double wn, int spec_index) const;

  //-----------------------------------------------------------------------
  /// The single scattering albedo for each layer, for the given wave number.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// state vector variables.
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \return Single scattering albedo for each layer. This is
  /// number_layer() in size
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> single_scattering_albedo_wrt_state_vector
  (double wn, int spec_index) const;

  //-----------------------------------------------------------------------
  /// The scattering moments for for each layer, for the given wave
  /// number.
  ///
  /// The scattering moments use the de Rooij convention for the 6
  /// scattering matrix element.
  ///
  /// The derivatives of the optical depth are with respect to the
  /// state vector variables.
  ///
  /// \param wn The wave number to calculate parameters for.
  /// \param spec_index The spectrometer index
  /// \param nummom Number of moments to include in
  ///           scatt_mom_each_layer, the default it to include all of
  ///           them.
  /// \param numscat Number of scattering matrix elements to include in
  ///           scatt_mom_each_layer, the default it to include all of
  ///           them.
  /// \return Scattering moments for each layer. This is
  ///         number_moment + 1 x number_layer() x number scattering
  ///         matrix elements
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 3> phase_function_moments_wrt_state_vector
  (double wn, int spec_index, int nummom = -1, int numscat = -1) const;

  //-----------------------------------------------------------------------
  /// The atmospheric thermal blackbody values per level.
  //-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> atmosphere_blackbody(double wn, int spec_index)
    const = 0;

  //-----------------------------------------------------------------------
  /// The surface thermal blackbody.
  //-----------------------------------------------------------------------
  
  virtual AutoDerivative<double> surface_blackbody(double wn, int spec_index)
    const = 0;

  //-----------------------------------------------------------------------
  /// This gives the optical properties value for a given wavenumber and
  /// spectral channel.
  //-----------------------------------------------------------------------
  
  virtual boost::shared_ptr<OpticalProperties> optical_properties
  (double wn, int spec_index) const = 0;
  
  //-----------------------------------------------------------------------
  /// This gives the values of the intermediate jacobian and that is used to
  /// translate the jacobians with respect to the RT parameters to a jacobian
  /// with respect to the state vector values.
  //-----------------------------------------------------------------------
  
  virtual blitz::Array<double, 3> intermediate_jacobian
  (double wn, int spec_index) const;

  //-----------------------------------------------------------------------
  /// Object that represents the ground surface. If null then there
  /// is no surface for this atmosphere.
  //-----------------------------------------------------------------------
  
  virtual const boost::shared_ptr<Ground> ground() const = 0;
  
  //-----------------------------------------------------------------------
  /// Reset timer
  //-----------------------------------------------------------------------
  
  virtual void reset_timer()
  {
    timer.reset_elapsed();
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(RtAtmosphere);
FP_EXPORT_OBSERVER_KEY(RtAtmosphere);
#endif
