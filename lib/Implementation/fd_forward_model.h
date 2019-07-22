#ifndef FD_FORWARD_MODEL_H
#define FD_FORWARD_MODEL_H
#include "forward_model.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This a forward model class that calcualtes jacobians using the
  finite difference method. This will approach is certainly
  much slower and is meant for testing analytic jacobains and
  for models that do not include analytic jacobians.

  Additionally this way of doing finitie difference is wasteful
  because some state vector elements might not require a full
  radiative transfer call.
*******************************************************************/

class FdForwardModel : public ForwardModel {
public:
  FdForwardModel(const boost::shared_ptr<ForwardModel>& Real_Forward_model,
  		 const boost::shared_ptr<StateVector>& Sv,
  		 const blitz::Array<double, 1>& Perturbation);
    virtual ~FdForwardModel() {}
    virtual boost::shared_ptr<StateVector> state_vector() const 
    { return statev; }
    virtual void print(std::ostream& Os) const { Os << "FdForwardModel"; }
    virtual Spectrum radiance(int Spec_index, bool Skip_jacobian = false) const;
    virtual int num_channels() const 
    {return real_fm->num_channels();}
    virtual SpectralDomain spectral_domain(int Spec_index) const
    {return real_fm->spectral_domain(Spec_index);}
private:
    boost::shared_ptr<ForwardModel> real_fm;
    boost::shared_ptr<StateVector> statev;
    blitz::Array<double, 1> perturb;
};
}
#endif
