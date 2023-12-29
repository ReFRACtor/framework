#ifndef FORWARD_MODEL_H
#define FORWARD_MODEL_H

#include "printable.h"
#include "spectrum.h"
#include "stacked_radiance_mixin.h"

namespace FullPhysics {
/****************************************************************//**
The forward model represents the encapsulation of modeling
spectra from an atmospheric state then applying instrument
specific effects to it.
*******************************************************************/
class ForwardModel : public StackedRadianceMixin {
public:
  virtual ~ForwardModel() {}

  //-----------------------------------------------------------------------
  /// This notifies the forward model that it should setup the grid
  //-----------------------------------------------------------------------

  virtual void setup_grid() = 0;

  //-----------------------------------------------------------------------
  /// The number of spectral channels associated with forward model.
  //-----------------------------------------------------------------------

  virtual int num_channels() const = 0;

  //-----------------------------------------------------------------------
  /// Spectral domain for the given spectral band. Note that this may be
  /// empty.
  //-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;

  //-----------------------------------------------------------------------
  /// Type preference for spectral domain. This may seem an odd thing to
  /// have a function for, but this is needed by ForwardModelOutput.
  //-----------------------------------------------------------------------

  virtual SpectralDomain::TypePreference spectral_domain_type_preference() const = 0;

  //-----------------------------------------------------------------------
  /// Spectrum for the given spectral band. Note that this may be empty.
  ///
  /// \param sensor_index Band to give value for
  /// \param skip_jacobian If true, don't do the Jacobian
  /// calculation. Often this is significantly faster to calculate.
  /// \return The set of radiances, possibly empty.
  //-----------------------------------------------------------------------

  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const = 0;

  virtual void print(std::ostream& Os) const { Os << "ForwardModel"; }

//-----------------------------------------------------------------------
/// We have some fairly nested object hierarchies. It can be useful to
/// be able to search this for things (e.g., which Pressure object is
/// used by a ForwardModel?). This returns a list of subobjects
/// "owned" by this object.
//-----------------------------------------------------------------------

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    return res;
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ForwardModel);
#endif
